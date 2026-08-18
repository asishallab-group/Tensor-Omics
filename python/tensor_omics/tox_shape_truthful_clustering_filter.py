r"""tox_shape_truthful_clustering_filter

# Shape Truthful Clustering (STC): Ensemble Filtering

`filter_ensembles`: decides which ensembles are eligible to be submitted to
`merge_to_super_ensembles` (see the sibling `tox_shape_truthful_clustering_reconciliation_impl`,
whose own `ensemble_reconciliation` is now a thin two-call orchestrator: this module first,
then that one). Composed of independent, individually testable per-criterion filters --
`filter_ensembles_by_stop_condition`, `filter_ensembles_by_dimension`,
`filter_ensembles_by_var_explained` -- each returning its own `eligible(n_ensembles)` mask
over the *same* ensembles, combined by `filter_ensembles` itself via a plain logical AND. A
criterion whose own threshold/allowed-set argument is omitted contributes an all-`True`
mask (no constraint from that criterion), so omitting every optional argument makes
`filter_ensembles` itself a true no-op (every ensemble eligible) -- see `misc/mod_STC.md`,
"Ensemble Reconciliation".

Filtering never alters `ensemble_identification`'s own output, nor does it remove an
ineligible ensemble from anywhere else this whole family reports it (points, the JSON's
`ensembles` array, CSV output, ...) -- only `merge_to_super_ensembles`'s own pairing/grouping
decision (and, downstream, `tox_stc_json`'s independently-computed `overlap_coefficient_matrix`,
which applies the identical mask for the same reason) ever sees an ineligible ensemble
excluded. No array copying/compaction anywhere in this module: every mask is exactly
`n_ensembles` long, over the same 1-indexed ensemble numbering everything else in this
family already uses.

Python binding, generated from tox_shape_truthful_clustering_filter. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.filter_ensembles_by_stop_condition_c.restype = None
_lib.filter_ensembles_by_stop_condition_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS')),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_ENSEMBLES_BY_STOP_CONDITION_ARGUMENTS = ("n_ensembles", "ensemble_stop_reason", "allowed_stop_reasons", "eligible", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_ENSEMBLES_BY_STOP_CONDITION_ARGUMENT_SOURCES = ("ensemble_stop_reason", None, None, None, None,)

_lib.filter_ensembles_by_dimension_c.restype = None
_lib.filter_ensembles_by_dimension_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_ENSEMBLES_BY_DIMENSION_ARGUMENTS = ("n_dimensions", "n_ensembles", "ensemble_d_final", "ensemble_has_final", "filter_dim_min", "filter_dim_max", "eligible", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_ENSEMBLES_BY_DIMENSION_ARGUMENT_SOURCES = (None, "ensemble_d_final", None, None, None, None, None, None,)

_lib.filter_ensembles_by_var_explained_c.restype = None
_lib.filter_ensembles_by_var_explained_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    nullable(ctypes.POINTER(ctypes.c_double)),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_ENSEMBLES_BY_VAR_EXPLAINED_ARGUMENTS = ("n_dimensions", "n_ensembles", "ensemble_S_final", "ensemble_d_final", "ensemble_k_final", "ensemble_has_final", "var_explained_min", "eligible", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_ENSEMBLES_BY_VAR_EXPLAINED_ARGUMENT_SOURCES = ("ensemble_S_final", "ensemble_S_final", None, None, None, None, None, None, None,)

_lib.filter_ensembles_c.restype = None
_lib.filter_ensembles_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS')),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_ENSEMBLES_ARGUMENTS = ("n_dimensions", "o", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_stop_reason", "allowed_stop_reasons", "filter_dim_min", "filter_dim_max", "var_explained_min", "eligible", "eligible_by_stop_condition", "eligible_by_dimension", "eligible_by_var_explained", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_ENSEMBLES_ARGUMENT_SOURCES = ("ensemble_U_history", "ensemble_U_history", "ensemble_U_history", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def filter_ensembles_by_stop_condition(
        ensemble_stop_reason,
        allowed_stop_reasons=None,
):
    r"""Ensemble eligibility by Stop Condition

    Parameters
    ----------
    ensemble_stop_reason : np.ndarray[np.int32] of shape (n_ensembles,)
        Per-ensemble Stop Condition, see `ensemble_identification`'s merged
        `ensemble_stop_reason` -- an index 1..4 into `allowed_stop_reasons` below, in the
        order `tox_shape_truthful_clustering_impl`'s own `STOP_REASON_MAX_SIZE` (1),
        `STOP_REASON_REJECTED_AFTER_STABLE` (2), `STOP_REASON_REJECTED_IMMEDIATELY` (3),
        `STOP_REASON_FIXED_POINT` (4) -- not imported by name here, to avoid a circular
        module dependency (the parent module already `use`s the reconciliation module,
        which `use`s this one)
        The minimum valid value is `1`.
        The maximum valid value is `4`.
    allowed_stop_reasons : np.ndarray[np.bool_] of shape (4,), optional
        Per-Stop-Condition eligibility, indexed as documented on `ensemble_stop_reason`
        above. Absent means no filtering (every Stop Condition allowed) -- deliberately
        nullable, not annotated with a generated default: the generator only evaluates
        constant *scalar* expressions for that annotation (`codegen_guide.md` section 5.5)

    Returns
    -------
    eligible : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
        Per-ensemble eligibility from this criterion alone
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_filter::filter_ensembles_by_stop_condition`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_stop_reason = np.ascontiguousarray(ensemble_stop_reason, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_stop_reason' must be an array of np.int32: {error}") from None
    if ensemble_stop_reason.ndim != 1:
        raise ValueError(f"'ensemble_stop_reason' must have 1 dimension, but has {ensemble_stop_reason.ndim}")
    if allowed_stop_reasons is not None:
        try:
            allowed_stop_reasons = np.ascontiguousarray(allowed_stop_reasons, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'allowed_stop_reasons' must be an array of np.bool_: {error}") from None
        if allowed_stop_reasons.ndim != 1:
            raise ValueError(f"'allowed_stop_reasons' must have 1 dimension, but has {allowed_stop_reasons.ndim}")

    # what the inputs already say, rather than asking for it again
    n_ensembles = ensemble_stop_reason.shape[0]

    # outputs and work arrays, which the caller never sees
    eligible = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.filter_ensembles_by_stop_condition_c(
        ctypes.byref(ctypes.c_int(n_ensembles)),
        ensemble_stop_reason,
        allowed_stop_reasons,
        eligible,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_ENSEMBLES_BY_STOP_CONDITION_ARGUMENTS, _FILTER_ENSEMBLES_BY_STOP_CONDITION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    eligible.flags.writeable = False

    return eligible

def filter_ensembles_by_dimension(
        n_dimensions,
        ensemble_d_final,
        ensemble_has_final,
        filter_dim_min=None,
        filter_dim_max=None,
):
    r"""Ensemble eligibility by final intrinsic dimension

    Parameters
    ----------
    n_dimensions : int
        Ambient dimension D
        The minimum valid value is `2`.
    ensemble_d_final : np.ndarray[np.int32] of shape (n_ensembles,)
        Each ensemble's final accepted intrinsic dimension, see
        `ensemble_final_observable`
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    ensemble_has_final : np.ndarray[np.bool_] of shape (n_ensembles,)
        Whether each ensemble has a final accepted state at all, see
        `ensemble_final_observable`
    filter_dim_min : int, optional
        Minimum tolerated final intrinsic dimension, inclusive
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    filter_dim_max : int, optional
        Maximum tolerated final intrinsic dimension, inclusive
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.

    Returns
    -------
    eligible : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
        Per-ensemble eligibility from this criterion alone
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_filter::filter_ensembles_by_dimension`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_d_final = np.ascontiguousarray(ensemble_d_final, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_d_final' must be an array of np.int32: {error}") from None
    if ensemble_d_final.ndim != 1:
        raise ValueError(f"'ensemble_d_final' must have 1 dimension, but has {ensemble_d_final.ndim}")
    try:
        ensemble_has_final = np.ascontiguousarray(ensemble_has_final, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_has_final' must be an array of np.bool_: {error}") from None
    if ensemble_has_final.ndim != 1:
        raise ValueError(f"'ensemble_has_final' must have 1 dimension, but has {ensemble_has_final.ndim}")

    # what the inputs already say, rather than asking for it again
    n_ensembles = ensemble_d_final.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if ensemble_has_final.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_has_final' has {ensemble_has_final.shape[0]} along axis 0, but "
            f"'ensemble_d_final' implies n_ensembles == {n_ensembles}"
        )

    # outputs and work arrays, which the caller never sees
    eligible = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.filter_ensembles_by_dimension_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        ensemble_d_final,
        ensemble_has_final,
        None if filter_dim_min is None else ctypes.byref(ctypes.c_int(filter_dim_min)),
        None if filter_dim_max is None else ctypes.byref(ctypes.c_int(filter_dim_max)),
        eligible,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_ENSEMBLES_BY_DIMENSION_ARGUMENTS, _FILTER_ENSEMBLES_BY_DIMENSION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    eligible.flags.writeable = False

    return eligible

def filter_ensembles_by_var_explained(
        ensemble_S_final,
        ensemble_d_final,
        ensemble_k_final,
        ensemble_has_final,
        var_explained_min=None,
):
    r"""Ensemble eligibility by final classical variance explained

    Parameters
    ----------
    ensemble_S_final : np.ndarray[np.float64] of shape (n_dimensions, n_ensembles,), column-major (order='F')
        Each ensemble's final accepted singular values, see `ensemble_final_observable`
    ensemble_d_final : np.ndarray[np.int32] of shape (n_ensembles,)
        Each ensemble's final accepted intrinsic dimension
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    ensemble_k_final : np.ndarray[np.int32] of shape (n_ensembles,)
        Each ensemble's final accepted size
        The minimum valid value is `0`.
    ensemble_has_final : np.ndarray[np.bool_] of shape (n_ensembles,)
        Whether each ensemble has a final accepted state at all
    var_explained_min : float, optional
        Minimum tolerated fraction of variance explained by the tangent subspace,
        inclusive
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    eligible : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
        Per-ensemble eligibility from this criterion alone
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_filter::filter_ensembles_by_var_explained`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_S_final = np.asfortranarray(ensemble_S_final, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_S_final' must be an array of np.float64: {error}") from None
    if ensemble_S_final.ndim != 2:
        raise ValueError(f"'ensemble_S_final' must have 2 dimensions, but has {ensemble_S_final.ndim}")
    try:
        ensemble_d_final = np.ascontiguousarray(ensemble_d_final, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_d_final' must be an array of np.int32: {error}") from None
    if ensemble_d_final.ndim != 1:
        raise ValueError(f"'ensemble_d_final' must have 1 dimension, but has {ensemble_d_final.ndim}")
    try:
        ensemble_k_final = np.ascontiguousarray(ensemble_k_final, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_k_final' must be an array of np.int32: {error}") from None
    if ensemble_k_final.ndim != 1:
        raise ValueError(f"'ensemble_k_final' must have 1 dimension, but has {ensemble_k_final.ndim}")
    try:
        ensemble_has_final = np.ascontiguousarray(ensemble_has_final, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_has_final' must be an array of np.bool_: {error}") from None
    if ensemble_has_final.ndim != 1:
        raise ValueError(f"'ensemble_has_final' must have 1 dimension, but has {ensemble_has_final.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = ensemble_S_final.shape[0]
    n_ensembles = ensemble_S_final.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if ensemble_d_final.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_d_final' has {ensemble_d_final.shape[0]} along axis 0, but "
            f"'ensemble_S_final' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_k_final.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_k_final' has {ensemble_k_final.shape[0]} along axis 0, but "
            f"'ensemble_S_final' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_has_final.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_has_final' has {ensemble_has_final.shape[0]} along axis 0, but "
            f"'ensemble_S_final' implies n_ensembles == {n_ensembles}"
        )

    # outputs and work arrays, which the caller never sees
    eligible = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.filter_ensembles_by_var_explained_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        ensemble_S_final,
        ensemble_d_final,
        ensemble_k_final,
        ensemble_has_final,
        None if var_explained_min is None else ctypes.byref(ctypes.c_double(var_explained_min)),
        eligible,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_ENSEMBLES_BY_VAR_EXPLAINED_ARGUMENTS, _FILTER_ENSEMBLES_BY_VAR_EXPLAINED_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    eligible.flags.writeable = False

    return eligible

def filter_ensembles(
        ensemble_U_history,
        ensemble_d_history,
        ensemble_S_history,
        ensemble_mu_history,
        ensemble_G_history,
        ensemble_k_history,
        ensemble_accepted_history,
        ensemble_stop_reason,
        allowed_stop_reasons=None,
        filter_dim_min=None,
        filter_dim_max=None,
        var_explained_min=None,
):
    r"""Combined ensemble eligibility for `merge_to_super_ensembles`

    Parameters
    ----------
    ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
        merged output
    ensemble_d_history : np.ndarray[np.int32] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing intrinsic dimensions
    ensemble_S_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing singular values
    ensemble_mu_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing centers
    ensemble_G_history : np.ndarray[np.float64] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing spectral gaps
    ensemble_k_history : np.ndarray[np.int32] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing sizes
    ensemble_accepted_history : np.ndarray[np.bool_] of shape (o, n_ensembles,), column-major (order='F')
        Whether the growth iteration retained in each history column was itself accepted
    ensemble_stop_reason : np.ndarray[np.int32] of shape (n_ensembles,)
        Per-ensemble Stop Condition, see `filter_ensembles_by_stop_condition_impl`
        The minimum valid value is `1`.
        The maximum valid value is `4`.
    allowed_stop_reasons : np.ndarray[np.bool_] of shape (4,), optional
        See `filter_ensembles_by_stop_condition_impl`
    filter_dim_min : int, optional
        See `filter_ensembles_by_dimension_impl`
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    filter_dim_max : int, optional
        See `filter_ensembles_by_dimension_impl`
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    var_explained_min : float, optional
        See `filter_ensembles_by_var_explained_impl`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        eligible : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            Combined eligibility: `True` only where all three per-criterion masks are
            A result is a value; call `.copy()` to obtain a modifiable array.
        eligible_by_stop_condition : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            See `filter_ensembles_by_stop_condition_impl`
            A result is a value; call `.copy()` to obtain a modifiable array.
        eligible_by_dimension : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            See `filter_ensembles_by_dimension_impl`
            A result is a value; call `.copy()` to obtain a modifiable array.
        eligible_by_var_explained : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            See `filter_ensembles_by_var_explained_impl`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_filter::filter_ensembles`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_U_history = np.asfortranarray(ensemble_U_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_U_history' must be an array of np.float64: {error}") from None
    if ensemble_U_history.ndim != 4:
        raise ValueError(f"'ensemble_U_history' must have 4 dimensions, but has {ensemble_U_history.ndim}")
    try:
        ensemble_d_history = np.asfortranarray(ensemble_d_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_d_history' must be an array of np.int32: {error}") from None
    if ensemble_d_history.ndim != 2:
        raise ValueError(f"'ensemble_d_history' must have 2 dimensions, but has {ensemble_d_history.ndim}")
    try:
        ensemble_S_history = np.asfortranarray(ensemble_S_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_S_history' must be an array of np.float64: {error}") from None
    if ensemble_S_history.ndim != 3:
        raise ValueError(f"'ensemble_S_history' must have 3 dimensions, but has {ensemble_S_history.ndim}")
    try:
        ensemble_mu_history = np.asfortranarray(ensemble_mu_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_mu_history' must be an array of np.float64: {error}") from None
    if ensemble_mu_history.ndim != 3:
        raise ValueError(f"'ensemble_mu_history' must have 3 dimensions, but has {ensemble_mu_history.ndim}")
    try:
        ensemble_G_history = np.asfortranarray(ensemble_G_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_G_history' must be an array of np.float64: {error}") from None
    if ensemble_G_history.ndim != 2:
        raise ValueError(f"'ensemble_G_history' must have 2 dimensions, but has {ensemble_G_history.ndim}")
    try:
        ensemble_k_history = np.asfortranarray(ensemble_k_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_k_history' must be an array of np.int32: {error}") from None
    if ensemble_k_history.ndim != 2:
        raise ValueError(f"'ensemble_k_history' must have 2 dimensions, but has {ensemble_k_history.ndim}")
    try:
        ensemble_accepted_history = np.asfortranarray(ensemble_accepted_history, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_accepted_history' must be an array of np.bool_: {error}") from None
    if ensemble_accepted_history.ndim != 2:
        raise ValueError(f"'ensemble_accepted_history' must have 2 dimensions, but has {ensemble_accepted_history.ndim}")
    try:
        ensemble_stop_reason = np.ascontiguousarray(ensemble_stop_reason, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_stop_reason' must be an array of np.int32: {error}") from None
    if ensemble_stop_reason.ndim != 1:
        raise ValueError(f"'ensemble_stop_reason' must have 1 dimension, but has {ensemble_stop_reason.ndim}")
    if allowed_stop_reasons is not None:
        try:
            allowed_stop_reasons = np.ascontiguousarray(allowed_stop_reasons, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'allowed_stop_reasons' must be an array of np.bool_: {error}") from None
        if allowed_stop_reasons.ndim != 1:
            raise ValueError(f"'allowed_stop_reasons' must have 1 dimension, but has {allowed_stop_reasons.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = ensemble_U_history.shape[0]
    o = ensemble_U_history.shape[2]
    n_ensembles = ensemble_U_history.shape[3]

    # Fortran cannot check that shared extents agree; this can
    if ensemble_S_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_mu_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_d_history.shape[0] != o:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_S_history.shape[1] != o:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_mu_history.shape[1] != o:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_G_history.shape[0] != o:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_k_history.shape[0] != o:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_accepted_history.shape[0] != o:
        raise ValueError(f"'ensemble_accepted_history' has {ensemble_accepted_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_d_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_S_history.shape[2] != n_ensembles:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[2]} along axis 2, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_mu_history.shape[2] != n_ensembles:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[2]} along axis 2, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_G_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_k_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_accepted_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_accepted_history' has {ensemble_accepted_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_stop_reason.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_stop_reason' has {ensemble_stop_reason.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies n_ensembles == {n_ensembles}"
        )

    # outputs and work arrays, which the caller never sees
    eligible = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_stop_condition = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_dimension = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_var_explained = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.filter_ensembles_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(o)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        ensemble_U_history,
        ensemble_d_history,
        ensemble_S_history,
        ensemble_mu_history,
        ensemble_G_history,
        ensemble_k_history,
        ensemble_accepted_history,
        ensemble_stop_reason,
        allowed_stop_reasons,
        None if filter_dim_min is None else ctypes.byref(ctypes.c_int(filter_dim_min)),
        None if filter_dim_max is None else ctypes.byref(ctypes.c_int(filter_dim_max)),
        None if var_explained_min is None else ctypes.byref(ctypes.c_double(var_explained_min)),
        eligible,
        eligible_by_stop_condition,
        eligible_by_dimension,
        eligible_by_var_explained,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FILTER_ENSEMBLES_ARGUMENTS, _FILTER_ENSEMBLES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    eligible.flags.writeable = False
    eligible_by_stop_condition.flags.writeable = False
    eligible_by_dimension.flags.writeable = False
    eligible_by_var_explained.flags.writeable = False

    return {
        "eligible": eligible,
        "eligible_by_stop_condition": eligible_by_stop_condition,
        "eligible_by_dimension": eligible_by_dimension,
        "eligible_by_var_explained": eligible_by_var_explained,
    }
