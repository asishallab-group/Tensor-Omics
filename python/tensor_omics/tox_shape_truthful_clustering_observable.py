r"""tox_shape_truthful_clustering_observable

# Shape Truthful Clustering (STC): Observable

`observable`: the tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble,
obtained from the economy-mode singular value decomposition (LAPACK `dgesdd`) of its
centered member vectors -- never an eigendecomposition of an explicitly formed
covariance matrix (see `misc/mod_STC.md`, "Numerical Linear Algebra"). `normal_error` and
`tangent_scales` are simple, dependency-free reductions over the eigenvalues `observable`
computes. See `misc/mod_STC.md`, SKG `observable`/`normal_error`/`tangent_scales`, for the
full algorithm definitions.

Python binding, generated from tox_shape_truthful_clustering_observable. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.normal_error_c.restype = None
_lib.normal_error_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMAL_ERROR_ARGUMENTS = ("d", "eigenvalues", "n_dimensions", "normal_error_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMAL_ERROR_ARGUMENT_SOURCES = (None, None, "eigenvalues", None, None,)

_lib.tangent_scales_c.restype = None
_lib.tangent_scales_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TANGENT_SCALES_ARGUMENTS = ("d", "eigenvalues", "n_dimensions", "tangent_scales_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_TANGENT_SCALES_ARGUMENT_SOURCES = ("tangent_scales_value", None, "eigenvalues", None, None,)

_lib.observable_c.restype = None
_lib.observable_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_OBSERVABLE_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "member_selection_mask", "n_selected_member", "U", "eigenvalues", "mu", "d", "G", "normal_error_value", "tangent_scales_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_OBSERVABLE_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, "member_selection_mask", None, None, None, None, None, None, None, None,)

_lib.ensemble_final_observable_c.restype = None
_lib.ensemble_final_observable_c.argtypes = (
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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ENSEMBLE_FINAL_OBSERVABLE_ARGUMENTS = ("n_dimensions", "o", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_U_final", "ensemble_d_final", "ensemble_S_final", "ensemble_mu_final", "ensemble_G_final", "ensemble_k_final", "ensemble_has_final", "ensemble_final_index", "ierr",)
#: For a derived argument, the one the caller passed it in
_ENSEMBLE_FINAL_OBSERVABLE_ARGUMENT_SOURCES = ("ensemble_U_history", "ensemble_U_history", "ensemble_U_history", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def normal_error(
        d,
        eigenvalues,
):
    r"""Mean squared residual of an ensemble's members off its tangent subspace

    Parameters
    ----------
    d : int
        Intrinsic (tangent) dimension of the ensemble
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,)
        Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        The minimum valid value is `0.0`.

    Returns
    -------
    normal_error_value : float
        Mean squared residual off the d-dimensional tangent subspace

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::normal_error`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        eigenvalues = np.ascontiguousarray(eigenvalues, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eigenvalues' must be an array of np.float64: {error}") from None
    if eigenvalues.ndim != 1:
        raise ValueError(f"'eigenvalues' must have 1 dimension, but has {eigenvalues.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = eigenvalues.shape[0]

    # outputs and work arrays, which the caller never sees
    normal_error_value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.normal_error_c(
        ctypes.byref(ctypes.c_int(d)),
        eigenvalues,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(normal_error_value),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMAL_ERROR_ARGUMENTS, _NORMAL_ERROR_ARGUMENT_SOURCES)

    return normal_error_value.value

def tangent_scales(
        d,
        eigenvalues,
):
    r"""Extent along each tangent direction of an ensemble's tangent subspace

    Parameters
    ----------
    d : int
        Intrinsic (tangent) dimension of the ensemble
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,)
        Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        The minimum valid value is `0.0`.

    Returns
    -------
    tangent_scales_value : np.ndarray[np.float64] of shape (d,), read-only
        Extent along each of the d tangent directions
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::tangent_scales`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        eigenvalues = np.ascontiguousarray(eigenvalues, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eigenvalues' must be an array of np.float64: {error}") from None
    if eigenvalues.ndim != 1:
        raise ValueError(f"'eigenvalues' must have 1 dimension, but has {eigenvalues.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = eigenvalues.shape[0]

    # outputs and work arrays, which the caller never sees
    tangent_scales_value = np.empty((d,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.tangent_scales_c(
        ctypes.byref(ctypes.c_int(d)),
        eigenvalues,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        tangent_scales_value,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TANGENT_SCALES_ARGUMENTS, _TANGENT_SCALES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    tangent_scales_value.flags.writeable = False

    return tangent_scales_value

def observable(
        vectors,
        member_selection_mask,
):
    r"""The tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    member_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Ensemble membership over the full dataset

    Returns
    -------
    dict
        with keys:

        U : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F'), read-only
            Tangent+normal basis, zero-padded beyond rank
            A result is a value; call `.copy()` to obtain a modifiable array.
        eigenvalues : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Covariance eigenvalues, descending, zero-padded beyond rank
            A result is a value; call `.copy()` to obtain a modifiable array.
        mu : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Ensemble center
            A result is a value; call `.copy()` to obtain a modifiable array.
        d : int
            Estimated intrinsic (tangent) dimension
        G : float
            Spectral gap at d
        normal_error_value : float
            Mean squared residual off the tangent subspace
        tangent_scales_value : np.ndarray[np.float64] of shape (n_dimensions,), read-only
            Extent along each tangent direction, zero-padded beyond d
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::observable`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vectors = np.asfortranarray(vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors' must be an array of np.float64: {error}") from None
    if vectors.ndim != 2:
        raise ValueError(f"'vectors' must have 2 dimensions, but has {vectors.ndim}")
    try:
        member_selection_mask = np.ascontiguousarray(member_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'member_selection_mask' must be an array of np.bool_: {error}") from None
    if member_selection_mask.ndim != 1:
        raise ValueError(f"'member_selection_mask' must have 1 dimension, but has {member_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_selected_member = int(member_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if member_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'member_selection_mask' has {member_selection_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    U = np.empty((n_dimensions, n_dimensions,), dtype=np.float64, order='F')
    eigenvalues = np.empty((n_dimensions,), dtype=np.float64, order='C')
    mu = np.empty((n_dimensions,), dtype=np.float64, order='C')
    d = ctypes.c_int(0)
    G = ctypes.c_double(0)
    normal_error_value = ctypes.c_double(0)
    tangent_scales_value = np.empty((n_dimensions,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.observable_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        member_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_member)),
        U,
        eigenvalues,
        mu,
        ctypes.byref(d),
        ctypes.byref(G),
        ctypes.byref(normal_error_value),
        tangent_scales_value,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _OBSERVABLE_ARGUMENTS, _OBSERVABLE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    U.flags.writeable = False
    eigenvalues.flags.writeable = False
    mu.flags.writeable = False
    tangent_scales_value.flags.writeable = False

    return {
        "U": U,
        "eigenvalues": eigenvalues,
        "mu": mu,
        "d": d.value,
        "G": G.value,
        "normal_error_value": normal_error_value.value,
        "tangent_scales_value": tangent_scales_value,
    }

def ensemble_final_observable(
        ensemble_U_history,
        ensemble_d_history,
        ensemble_S_history,
        ensemble_mu_history,
        ensemble_G_history,
        ensemble_k_history,
        ensemble_accepted_history,
):
    r"""Each ensemble's final *accepted* growth-history state

    Parameters
    ----------
    ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
        merged `ensemble_U_history`
    ensemble_d_history : np.ndarray[np.int32] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing intrinsic dimensions
    ensemble_S_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing singular values
    ensemble_mu_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing centers
    ensemble_G_history : np.ndarray[np.float64] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing spectral gaps
    ensemble_k_history : np.ndarray[np.int32] of shape (o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing sizes; 0 marks an unpopulated column
    ensemble_accepted_history : np.ndarray[np.bool_] of shape (o, n_ensembles,), column-major (order='F')
        Whether the growth iteration retained in each history column was itself
        accepted -- see this kernel's own summary above

    Returns
    -------
    dict
        with keys:

        ensemble_U_final : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, n_ensembles,), column-major (order='F'), read-only
            Each ensemble's final accepted tangent+normal basis; zero when
            `ensemble_has_final` is `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_d_final : np.ndarray[np.int32] of shape (n_ensembles,), read-only
            Each ensemble's final accepted intrinsic dimension; zero when
            `ensemble_has_final` is `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_S_final : np.ndarray[np.float64] of shape (n_dimensions, n_ensembles,), column-major (order='F'), read-only
            Each ensemble's final accepted singular values; zero when
            `ensemble_has_final` is `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_mu_final : np.ndarray[np.float64] of shape (n_dimensions, n_ensembles,), column-major (order='F'), read-only
            Each ensemble's final accepted center; zero when `ensemble_has_final` is
            `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_G_final : np.ndarray[np.float64] of shape (n_ensembles,), read-only
            Each ensemble's final accepted spectral gap; zero when `ensemble_has_final` is
            `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_k_final : np.ndarray[np.int32] of shape (n_ensembles,), read-only
            Each ensemble's final accepted size; zero when `ensemble_has_final` is
            `False` for that ensemble
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_has_final : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            Whether any history column at all qualifies as this ensemble's final accepted
            state -- see this kernel's own summary above for the (rare) `False` cases
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_final_index : np.ndarray[np.int32] of shape (n_ensembles,), read-only
            The history column each `_final` output was sliced from (0 when
            `ensemble_has_final` is `False`) -- also, since every column 1..this index is
            itself guaranteed accepted (only ever the single *last* populated column can be
            the rejected candidate this kernel's own summary describes), this doubles as the
            count of genuinely accepted, plottable history columns, for callers (e.g.
            `tox_stc_json`'s own `observable_history`) that need to iterate the whole
            trailing window, not just its final entry
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable::ensemble_final_observable`, whose argument names are
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

    # outputs and work arrays, which the caller never sees
    ensemble_U_final = np.empty((n_dimensions, n_dimensions, n_ensembles,), dtype=np.float64, order='F')
    ensemble_d_final = np.empty((n_ensembles,), dtype=np.int32, order='C')
    ensemble_S_final = np.empty((n_dimensions, n_ensembles,), dtype=np.float64, order='F')
    ensemble_mu_final = np.empty((n_dimensions, n_ensembles,), dtype=np.float64, order='F')
    ensemble_G_final = np.empty((n_ensembles,), dtype=np.float64, order='C')
    ensemble_k_final = np.empty((n_ensembles,), dtype=np.int32, order='C')
    ensemble_has_final = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ensemble_final_index = np.empty((n_ensembles,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.ensemble_final_observable_c(
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
        ensemble_U_final,
        ensemble_d_final,
        ensemble_S_final,
        ensemble_mu_final,
        ensemble_G_final,
        ensemble_k_final,
        ensemble_has_final,
        ensemble_final_index,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ENSEMBLE_FINAL_OBSERVABLE_ARGUMENTS, _ENSEMBLE_FINAL_OBSERVABLE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ensemble_U_final.flags.writeable = False
    ensemble_d_final.flags.writeable = False
    ensemble_S_final.flags.writeable = False
    ensemble_mu_final.flags.writeable = False
    ensemble_G_final.flags.writeable = False
    ensemble_k_final.flags.writeable = False
    ensemble_has_final.flags.writeable = False
    ensemble_final_index.flags.writeable = False

    return {
        "ensemble_U_final": ensemble_U_final,
        "ensemble_d_final": ensemble_d_final,
        "ensemble_S_final": ensemble_S_final,
        "ensemble_mu_final": ensemble_mu_final,
        "ensemble_G_final": ensemble_G_final,
        "ensemble_k_final": ensemble_k_final,
        "ensemble_has_final": ensemble_has_final,
        "ensemble_final_index": ensemble_final_index,
    }
