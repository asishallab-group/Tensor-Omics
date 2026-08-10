"""tox_shape_truthful_clustering

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shape_truthful_clustering. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.ensemble_identification_c.restype = None
_lib.ensemble_identification_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ENSEMBLE_IDENTIFICATION_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_index", "k_min", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "o", "final_ensemble_mask", "stop_reason", "growth_radius", "U_history", "S_history", "d_history", "G_history", "mu_history", "k_history", "accepted_history", "member_added_at_step", "low_confidence_mask", "U_first", "d_first", "ierr",)
#: For a derived argument, the one the caller passed it in
_ENSEMBLE_IDENTIFICATION_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, None, None, None, None, None, None, None, "U_history", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.ensemble_identification_merged_c.restype = None
_lib.ensemble_identification_merged_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ENSEMBLE_IDENTIFICATION_MERGED_ARGUMENTS = ("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_selection_mask", "n_selected_seed", "k_min", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "o", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_member_added_at_step", "ensemble_low_confidence_masks", "ensemble_U_first", "ensemble_d_first", "ierr",)
#: For a derived argument, the one the caller passed it in
_ENSEMBLE_IDENTIFICATION_MERGED_ARGUMENT_SOURCES = (None, "vectors", "vectors", None, None, None, "ensemble_masks", None, None, None, None, None, None, None, "ensemble_U_history", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def ensemble_identification(
        vectors,
        kd_indices,
        dimension_order,
        seed_index,
        chordal_dist_max_as_prcnt_of_range,
        d_max,
        G_max,
        RMSE_change_max,
        o,
        k_min=30,
        f_max=0.95,
        a=2,
):
    r"""Grow and track a single ensemble from one seed until a Stop Condition is reached

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        Pre-built k-d tree index over `vectors`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    seed_index : int
        Index into `vectors`/`kd_indices` of the seed to grow an ensemble around
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    k_min : int, optional, default 30
        Neighborhood size for this seed's growth radius, see `calc_ensemble_growth_radius`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
        The default value is `30`.
    chordal_dist_max_as_prcnt_of_range : float
        Maximum tolerated chordal distance between tangent bases, as a fraction of its
        own [0, sqrt(d)] range, see `accept_ensemble`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    d_max : int
        Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
        The minimum valid value is `0`.
    G_max : float
        Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
        The minimum valid value is `0.0`.
    RMSE_change_max : float
        Maximum tolerated |log(RMSE_tp1/RMSE_t)|, see `accept_ensemble`
        The minimum valid value is `0.0`.
    f_max : float, optional, default 0.95
        Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
        The minimum valid value is `above(0.0)`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.
    a : int, optional, default 2
        Minimum accepted-iteration count for a later rejection to count as "stable", see
        Stop Condition 2
        The minimum valid value is `1`.
        The default value is `2`.
    o : int
        Trailing observable-history window depth (`misc/mod_STC.md` suggests 10 as a
        sensible default). Always required, never optional with an auto-applied
        default here: a Fortran array bound cannot depend on a possibly-absent
        optional dummy, and this argument sizes every history output below.
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        final_ensemble_mask : np.ndarray[np.bool_] of shape (n_vectors,), read-only
            The last accepted ensemble's membership. All `False` when `stop_reason` is
            `STOP_REASON_MAX_SIZE` -- see Stop Condition 1.
            A result is a value; call `.copy()` to obtain a modifiable array.
        stop_reason : int
            Which Stop Condition ended growth: one of
            ``STOP_REASON_MAX_SIZE``,
            ``STOP_REASON_REJECTED_AFTER_STABLE``,
            ``STOP_REASON_REJECTED_IMMEDIATELY``, or
            ``STOP_REASON_FIXED_POINT`` --
            or ``STOP_REASON_ERROR`` if
            `ierr` is non-zero, in which case every other output for this seed is undefined.
        growth_radius : float
            This seed's growth radius, see `calc_ensemble_growth_radius`
        U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o,), column-major (order='F'), read-only
            Trailing tangent+normal bases, one per retained iteration, oldest to newest;
            zero beyond the number of iterations actually retained, see `k_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        S_history : np.ndarray[np.float64] of shape (n_dimensions, o,), column-major (order='F'), read-only
            Trailing singular values -- not eigenvalues, see "Output" in `misc/mod_STC.md`
            -- zero-padded beyond rank and beyond the number of retained iterations
            A result is a value; call `.copy()` to obtain a modifiable array.
        d_history : np.ndarray[np.int32] of shape (o,), read-only
            Trailing intrinsic dimensions, one per retained iteration
            A result is a value; call `.copy()` to obtain a modifiable array.
        G_history : np.ndarray[np.float64] of shape (o,), read-only
            Trailing spectral gaps, one per retained iteration
            A result is a value; call `.copy()` to obtain a modifiable array.
        mu_history : np.ndarray[np.float64] of shape (n_dimensions, o,), column-major (order='F'), read-only
            Trailing ensemble centers, one per retained iteration
            A result is a value; call `.copy()` to obtain a modifiable array.
        k_history : np.ndarray[np.int32] of shape (o,), read-only
            Trailing ensemble sizes, one per retained iteration. 0 marks a column beyond
            the number of iterations actually retained -- a real ensemble size is always
            at least 1.
            A result is a value; call `.copy()` to obtain a modifiable array.
        accepted_history : np.ndarray[np.bool_] of shape (o,), read-only
            Whether the growth iteration retained in the corresponding column was
            accepted. Iteration 1 (the bootstrap step) is always `True` by convention.
            The single most recent column is `False` when, and only when, growth
            stopped via `STOP_REASON_REJECTED_AFTER_STABLE` or
            `STOP_REASON_REJECTED_IMMEDIATELY` -- see the module-level note above.
            A result is a value; call `.copy()` to obtain a modifiable array.
        member_added_at_step : np.ndarray[np.int32] of shape (n_vectors,), read-only
            `MEMBER_ADDED_AT_STEP_NON_MEMBER` for non-members, `MEMBER_ADDED_AT_STEP_SEED`
            for the seed itself, the growth-iteration index at which each other member
            joined otherwise
            A result is a value; call `.copy()` to obtain a modifiable array.
        low_confidence_mask : np.ndarray[np.bool_] of shape (n_vectors,), read-only
            Membership from this seed's iteration 1 (the unconditional bootstrap
            grow_ensemble+observable call), reported regardless of stop_reason -- including
            when stop_reason is STOP_REASON_MAX_SIZE, for which final_ensemble_mask is
            all-False All-False here too whenever iteration 1 itself never produced a
            genuine observable (an isolated seed, or a seed whose very first growth step
            already exceeds f_max*N) -- see "Ensemble identification", "Output" in
            misc/mod_STC.md
            A result is a value; call `.copy()` to obtain a modifiable array.
        U_first : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F'), read-only
            Tangent+normal basis at the bootstrap iteration (iteration 1), retained for the
            whole growth, never evicted by the trailing o-window above -- see
            `accept_ensemble`'s tangent-space-drift criterion in misc/mod_STC.md. All-zero
            whenever iteration 1 itself never produced a genuine observable, same condition
            as `low_confidence_mask` above.
            A result is a value; call `.copy()` to obtain a modifiable array.
        d_first : int
            Intrinsic dimension at the bootstrap iteration, see `U_first`. Zero under the
            same all-zero condition as `U_first`.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering::ensemble_identification`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_vectors:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    final_ensemble_mask = np.empty((n_vectors,), dtype=np.bool_, order='C')
    stop_reason = ctypes.c_int(0)
    growth_radius = ctypes.c_double(0)
    U_history = np.empty((n_dimensions, n_dimensions, o,), dtype=np.float64, order='F')
    S_history = np.empty((n_dimensions, o,), dtype=np.float64, order='F')
    d_history = np.empty((o,), dtype=np.int32, order='C')
    G_history = np.empty((o,), dtype=np.float64, order='C')
    mu_history = np.empty((n_dimensions, o,), dtype=np.float64, order='F')
    k_history = np.empty((o,), dtype=np.int32, order='C')
    accepted_history = np.empty((o,), dtype=np.bool_, order='C')
    member_added_at_step = np.empty((n_vectors,), dtype=np.int32, order='C')
    low_confidence_mask = np.empty((n_vectors,), dtype=np.bool_, order='C')
    U_first = np.empty((n_dimensions, n_dimensions,), dtype=np.float64, order='F')
    d_first = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.ensemble_identification_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        ctypes.byref(ctypes.c_int(seed_index)),
        ctypes.byref(ctypes.c_int(k_min)),
        ctypes.byref(ctypes.c_double(chordal_dist_max_as_prcnt_of_range)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(ctypes.c_double(RMSE_change_max)),
        ctypes.byref(ctypes.c_double(f_max)),
        ctypes.byref(ctypes.c_int(a)),
        ctypes.byref(ctypes.c_int(o)),
        final_ensemble_mask,
        ctypes.byref(stop_reason),
        ctypes.byref(growth_radius),
        U_history,
        S_history,
        d_history,
        G_history,
        mu_history,
        k_history,
        accepted_history,
        member_added_at_step,
        low_confidence_mask,
        U_first,
        ctypes.byref(d_first),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ENSEMBLE_IDENTIFICATION_ARGUMENTS, _ENSEMBLE_IDENTIFICATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    final_ensemble_mask.flags.writeable = False
    U_history.flags.writeable = False
    S_history.flags.writeable = False
    d_history.flags.writeable = False
    G_history.flags.writeable = False
    mu_history.flags.writeable = False
    k_history.flags.writeable = False
    accepted_history.flags.writeable = False
    member_added_at_step.flags.writeable = False
    low_confidence_mask.flags.writeable = False
    U_first.flags.writeable = False

    return {
        "final_ensemble_mask": final_ensemble_mask,
        "stop_reason": stop_reason.value,
        "growth_radius": growth_radius.value,
        "U_history": U_history,
        "S_history": S_history,
        "d_history": d_history,
        "G_history": G_history,
        "mu_history": mu_history,
        "k_history": k_history,
        "accepted_history": accepted_history,
        "member_added_at_step": member_added_at_step,
        "low_confidence_mask": low_confidence_mask,
        "U_first": U_first,
        "d_first": d_first.value,
    }

def ensemble_identification_merged(
        vectors,
        kd_indices,
        dimension_order,
        seed_selection_mask,
        chordal_dist_max_as_prcnt_of_range,
        d_max,
        G_max,
        RMSE_change_max,
        o,
        k_min=30,
        f_max=0.95,
        a=2,
):
    r"""Run ensemble_identification once per seed and assemble the merged, per-ensemble output arrays

    Parameters
    ----------
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    kd_indices : np.ndarray[np.int32] of shape (n_vectors,)
        Pre-built k-d tree index over `vectors`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors`.
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order used to build `kd_indices`
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.
    seed_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Seed selection, see `seeds`
    k_min : int, optional, default 30
        Neighborhood size for each seed's growth radius, see `calc_ensemble_growth_radius`
        The minimum valid value is `1`.
        The maximum valid value is `n_vectors - 1`.
        The default value is `30`.
    chordal_dist_max_as_prcnt_of_range : float
        Maximum tolerated chordal distance between tangent bases, as a fraction of its
        own [0, sqrt(d)] range, see `accept_ensemble`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    d_max : int
        Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
        The minimum valid value is `0`.
    G_max : float
        Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
        The minimum valid value is `0.0`.
    RMSE_change_max : float
        Maximum tolerated |log(RMSE_tp1/RMSE_t)|, see `accept_ensemble`
        The minimum valid value is `0.0`.
    f_max : float, optional, default 0.95
        Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
        The minimum valid value is `above(0.0)`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.
    a : int, optional, default 2
        Minimum accepted-iteration count for a later rejection to count as "stable", see
        Stop Condition 2
        The minimum valid value is `1`.
        The default value is `2`.
    o : int
        Trailing observable-history window depth, see `ensemble_identification`. Always
        required, for the same reason as there: it sizes every history output below.
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble accepted membership, one column per seed, see `final_ensemble_mask`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_stop_reason : np.ndarray[np.int32] of shape (n_selected_seed,), read-only
            Per-ensemble Stop Condition, see `ensemble_identification`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_growth_radii : np.ndarray[np.float64] of shape (n_selected_seed,), read-only
            Per-ensemble growth radius, see "Local Radius Identification"
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing tangent+normal bases, see `U_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_S_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing singular values, see `S_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_d_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing intrinsic dimensions, see `d_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_G_history : np.ndarray[np.float64] of shape (o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing spectral gaps, see `G_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_mu_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing centers, see `mu_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_k_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing sizes, see `k_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_accepted_history : np.ndarray[np.bool_] of shape (o, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble trailing accepted flags, see `accepted_history`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_member_added_at_step : np.ndarray[np.int32] of shape (n_vectors, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble growth-iteration-joined bookkeeping, see `member_added_at_step`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_low_confidence_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble iteration-1 fallback membership, see `low_confidence_mask`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_U_first : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, n_selected_seed,), column-major (order='F'), read-only
            Per-ensemble bootstrap-iteration tangent+normal basis, see `U_first`
            A result is a value; call `.copy()` to obtain a modifiable array.
        ensemble_d_first : np.ndarray[np.int32] of shape (n_selected_seed,), read-only
            Per-ensemble bootstrap-iteration intrinsic dimension, see `d_first`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering::ensemble_identification_merged`, whose argument names are
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
        kd_indices = np.ascontiguousarray(kd_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'kd_indices' must be an array of np.int32: {error}") from None
    if kd_indices.ndim != 1:
        raise ValueError(f"'kd_indices' must have 1 dimension, but has {kd_indices.ndim}")
    try:
        dimension_order = np.ascontiguousarray(dimension_order, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dimension_order' must be an array of np.int32: {error}") from None
    if dimension_order.ndim != 1:
        raise ValueError(f"'dimension_order' must have 1 dimension, but has {dimension_order.ndim}")
    try:
        seed_selection_mask = np.ascontiguousarray(seed_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'seed_selection_mask' must be an array of np.bool_: {error}") from None
    if seed_selection_mask.ndim != 1:
        raise ValueError(f"'seed_selection_mask' must have 1 dimension, but has {seed_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_selected_seed = int(seed_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if dimension_order.shape[0] != n_dimensions:
        raise ValueError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if kd_indices.shape[0] != n_vectors:
        raise ValueError(f"'kd_indices' has {kd_indices.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if seed_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'seed_selection_mask' has {seed_selection_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    ensemble_masks = np.empty((n_vectors, n_selected_seed,), dtype=np.bool_, order='F')
    ensemble_stop_reason = np.empty((n_selected_seed,), dtype=np.int32, order='C')
    ensemble_growth_radii = np.empty((n_selected_seed,), dtype=np.float64, order='C')
    ensemble_U_history = np.empty((n_dimensions, n_dimensions, o, n_selected_seed,), dtype=np.float64, order='F')
    ensemble_S_history = np.empty((n_dimensions, o, n_selected_seed,), dtype=np.float64, order='F')
    ensemble_d_history = np.empty((o, n_selected_seed,), dtype=np.int32, order='F')
    ensemble_G_history = np.empty((o, n_selected_seed,), dtype=np.float64, order='F')
    ensemble_mu_history = np.empty((n_dimensions, o, n_selected_seed,), dtype=np.float64, order='F')
    ensemble_k_history = np.empty((o, n_selected_seed,), dtype=np.int32, order='F')
    ensemble_accepted_history = np.empty((o, n_selected_seed,), dtype=np.bool_, order='F')
    ensemble_member_added_at_step = np.empty((n_vectors, n_selected_seed,), dtype=np.int32, order='F')
    ensemble_low_confidence_masks = np.empty((n_vectors, n_selected_seed,), dtype=np.bool_, order='F')
    ensemble_U_first = np.empty((n_dimensions, n_dimensions, n_selected_seed,), dtype=np.float64, order='F')
    ensemble_d_first = np.empty((n_selected_seed,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.ensemble_identification_merged_c(
        vectors,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        kd_indices,
        dimension_order,
        seed_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_seed)),
        ctypes.byref(ctypes.c_int(k_min)),
        ctypes.byref(ctypes.c_double(chordal_dist_max_as_prcnt_of_range)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(ctypes.c_double(RMSE_change_max)),
        ctypes.byref(ctypes.c_double(f_max)),
        ctypes.byref(ctypes.c_int(a)),
        ctypes.byref(ctypes.c_int(o)),
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_growth_radii,
        ensemble_U_history,
        ensemble_S_history,
        ensemble_d_history,
        ensemble_G_history,
        ensemble_mu_history,
        ensemble_k_history,
        ensemble_accepted_history,
        ensemble_member_added_at_step,
        ensemble_low_confidence_masks,
        ensemble_U_first,
        ensemble_d_first,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ENSEMBLE_IDENTIFICATION_MERGED_ARGUMENTS, _ENSEMBLE_IDENTIFICATION_MERGED_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ensemble_masks.flags.writeable = False
    ensemble_stop_reason.flags.writeable = False
    ensemble_growth_radii.flags.writeable = False
    ensemble_U_history.flags.writeable = False
    ensemble_S_history.flags.writeable = False
    ensemble_d_history.flags.writeable = False
    ensemble_G_history.flags.writeable = False
    ensemble_mu_history.flags.writeable = False
    ensemble_k_history.flags.writeable = False
    ensemble_accepted_history.flags.writeable = False
    ensemble_member_added_at_step.flags.writeable = False
    ensemble_low_confidence_masks.flags.writeable = False
    ensemble_U_first.flags.writeable = False
    ensemble_d_first.flags.writeable = False

    return {
        "ensemble_masks": ensemble_masks,
        "ensemble_stop_reason": ensemble_stop_reason,
        "ensemble_growth_radii": ensemble_growth_radii,
        "ensemble_U_history": ensemble_U_history,
        "ensemble_S_history": ensemble_S_history,
        "ensemble_d_history": ensemble_d_history,
        "ensemble_G_history": ensemble_G_history,
        "ensemble_mu_history": ensemble_mu_history,
        "ensemble_k_history": ensemble_k_history,
        "ensemble_accepted_history": ensemble_accepted_history,
        "ensemble_member_added_at_step": ensemble_member_added_at_step,
        "ensemble_low_confidence_masks": ensemble_low_confidence_masks,
        "ensemble_U_first": ensemble_U_first,
        "ensemble_d_first": ensemble_d_first,
    }
