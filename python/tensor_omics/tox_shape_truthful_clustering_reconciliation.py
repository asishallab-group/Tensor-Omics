r"""tox_shape_truthful_clustering_reconciliation

# Shape Truthful Clustering (STC): Ensemble Reconciliation

`ensemble_reconciliation`: identifies intersecting ensembles from Ensemble Identification's
merged `ensemble_masks` output and, depending on `mode`, either just reports intersecting
pairs or groups transitively-intersecting ensembles into "super-ensembles" via a union-find
over the pairwise intersection graph. See `misc/mod_STC.md`, "Ensemble Reconciliation", for
the full algorithm definition. Does not alter Ensemble Identification's own result -- this
module only reports and groups, on the side.

A thin, two-call orchestrator over its own two sibling kernels: first
:func:`tensor_omics.filter_ensembles` (which
ensembles are even eligible to contribute a pair, by Stop Condition/final dimension/final
variance explained), then this module's own `merge_to_super_ensembles_impl` (the actual
pairwise-intersection/union-find grouping, over eligible ensembles only). Splitting these
into two independently testable, independently reusable kernels -- rather than one kernel
that both decides eligibility and merges -- is a deliberate design choice: eligibility is a
statement about *individual* ensembles (their own Stop Condition/geometry), merging is a
statement about *pairs*, and conflating the two made every new filtering criterion require
touching the same monolithic merge logic. See `misc/mod_STC.md`'s own rationale for the
split.

Python binding, generated from tox_shape_truthful_clustering_reconciliation. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.ensemble_reconciliation_c.restype = None
_lib.ensemble_reconciliation_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
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
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS')),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ENSEMBLE_RECONCILIATION_ARGUMENTS = ("ensemble_masks", "ensemble_stop_reason", "n_dimensions", "n_vectors", "n_ensembles", "ensemble_U_history", "ensemble_d_history", "ensemble_S_history", "ensemble_mu_history", "ensemble_G_history", "ensemble_k_history", "ensemble_accepted_history", "o", "mode", "min_overlap_coefficient", "report_overlap_coefficient", "allowed_stop_reasons", "filter_dim_min", "filter_dim_max", "var_explained_min", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_overlap_coefficient", "eligible", "eligible_by_stop_condition", "eligible_by_dimension", "eligible_by_var_explained", "ierr",)
#: For a derived argument, the one the caller passed it in
_ENSEMBLE_RECONCILIATION_ARGUMENT_SOURCES = (None, None, "ensemble_U_history", "ensemble_masks", "ensemble_masks", None, None, None, None, None, None, None, "ensemble_U_history", None, None, None, None, None, None, None, "super_ensembles", None, None, None, None, None, None, None, None,)

_lib.merge_to_super_ensembles_c.restype = None
_lib.merge_to_super_ensembles_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MERGE_TO_SUPER_ENSEMBLES_ARGUMENTS = ("ensemble_masks", "eligible", "n_vectors", "n_ensembles", "mode", "min_overlap_coefficient", "report_overlap_coefficient", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_overlap_coefficient", "ierr",)
#: For a derived argument, the one the caller passed it in
_MERGE_TO_SUPER_ENSEMBLES_ARGUMENT_SOURCES = (None, None, "ensemble_masks", "ensemble_masks", None, None, None, "super_ensembles", None, None, None, None,)

def ensemble_reconciliation(
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_U_history,
        ensemble_d_history,
        ensemble_S_history,
        ensemble_mu_history,
        ensemble_G_history,
        ensemble_k_history,
        ensemble_accepted_history,
        max_group_size,
        mode='report',
        min_overlap_coefficient=0.9,
        report_overlap_coefficient=False,
        allowed_stop_reasons=None,
        filter_dim_min=None,
        filter_dim_max=None,
        var_explained_min=None,
):
    r"""Filter eligible ensembles, then group/report their intersections

    Parameters
    ----------
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_ensembles,), column-major (order='F')
        Per-ensemble membership, see Ensemble Identification's merged output
    ensemble_stop_reason : np.ndarray[np.int32] of shape (n_ensembles,)
        Per-ensemble Stop Condition, see `filter_ensembles_impl`
        The minimum valid value is `1`.
        The maximum valid value is `4`.
    ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_ensembles,), column-major (order='F')
        Per-ensemble trailing tangent+normal bases, see Ensemble Identification's merged
        output
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
    mode : str, one of 'report' | 'merge_overlap_coefficient' | 'merge_any', optional, default 'report'
        How intersections are processed

        The default value is `'report'`.
    min_overlap_coefficient : float, optional, default 0.9
        Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
        \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
        ``'merge_overlap_coefficient'``;
        ignored in every other mode
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    report_overlap_coefficient : bool, optional, default False
        Whether to compute and return `super_ensembles_overlap_coefficient` at all --
        see `merge_to_super_ensembles_impl`'s own note on this being guarded, not
        unconditional
        The default value is `False`.
    allowed_stop_reasons : np.ndarray[np.bool_] of shape (4,), optional
        See `tox_shape_truthful_clustering_filter_impl`'s own
        `filter_ensembles_by_stop_condition_impl`
    filter_dim_min : int, optional
        See `tox_shape_truthful_clustering_filter_impl`'s own
        `filter_ensembles_by_dimension_impl`
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    filter_dim_max : int, optional
        See `tox_shape_truthful_clustering_filter_impl`'s own
        `filter_ensembles_by_dimension_impl`
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    var_explained_min : float, optional
        See `tox_shape_truthful_clustering_filter_impl`'s own
        `filter_ensembles_by_var_explained_impl`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    max_group_size : int
        Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
        can hold; sizes its row dimension. `misc/mod_STC.md` suggests
        $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
        optional with an auto-applied default here, for the same reason as
        `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
        possibly-absent optional dummy, and a runtime-dependent value like
        $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
        default would need to be either.
        The minimum valid value is `2`.
        The maximum valid value is `n_ensembles`.

    Returns
    -------
    dict
        with keys:

        super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_ensembles*(n_ensembles-1),), column-major (order='F'), read-only
            One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            the group's actual size, and 0 in every row of an unused trailing column beyond
            `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            ``'report'``'s
            own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            intersects) -- deliberately not divided by 2: the generator translates this
            specification expression close to verbatim into the Python/R bindings, where
            `/` on two integers is true division, not Fortran's own truncating integer
            division, so a literal `/2` here breaks the generated Python binding (a `float`
            where `np.empty`'s shape wants an `int`); see `misc/code_gen_footgun.md`. A
            safe, if looser, upper bound for modes 2 and 3 too, whose groups can never
            outnumber mode 1's own worst case.
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_super_ensembles : int
            Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            actually filled
        super_ensembles_overlap_coefficient : np.ndarray[np.float64] of shape (max_group_size-1, n_ensembles*(n_ensembles-1),), column-major (order='F'), read-only
            Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            `report_overlap_coefficient` was requested -- see the note above.
            A result is a value; call `.copy()` to obtain a modifiable array.
        eligible : np.ndarray[np.bool_] of shape (n_ensembles,), read-only
            Combined per-ensemble eligibility actually used for merging above -- see
            `filter_ensembles_impl`. Ineligible ensembles are otherwise untouched: they
            are never removed from `ensemble_masks` or anything else this whole family
            reports, only excluded from contributing a pair here.
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
    Generated from the Fortran procedure `tox_shape_truthful_clustering_reconciliation::ensemble_reconciliation`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_masks = np.asfortranarray(ensemble_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_masks' must be an array of np.bool_: {error}") from None
    if ensemble_masks.ndim != 2:
        raise ValueError(f"'ensemble_masks' must have 2 dimensions, but has {ensemble_masks.ndim}")
    try:
        ensemble_stop_reason = np.ascontiguousarray(ensemble_stop_reason, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_stop_reason' must be an array of np.int32: {error}") from None
    if ensemble_stop_reason.ndim != 1:
        raise ValueError(f"'ensemble_stop_reason' must have 1 dimension, but has {ensemble_stop_reason.ndim}")
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
    mode = np.array([str(mode).lower().encode().ljust(25)], dtype="S25")
    if allowed_stop_reasons is not None:
        try:
            allowed_stop_reasons = np.ascontiguousarray(allowed_stop_reasons, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'allowed_stop_reasons' must be an array of np.bool_: {error}") from None
        if allowed_stop_reasons.ndim != 1:
            raise ValueError(f"'allowed_stop_reasons' must have 1 dimension, but has {allowed_stop_reasons.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = ensemble_U_history.shape[0]
    n_vectors = ensemble_masks.shape[0]
    n_ensembles = ensemble_masks.shape[1]
    o = ensemble_U_history.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if ensemble_S_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_mu_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_stop_reason.shape[0] != n_ensembles:
        raise ValueError(f"'ensemble_stop_reason' has {ensemble_stop_reason.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_U_history.shape[3] != n_ensembles:
        raise ValueError(f"'ensemble_U_history' has {ensemble_U_history.shape[3]} along axis 3, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_d_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_S_history.shape[2] != n_ensembles:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_mu_history.shape[2] != n_ensembles:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_G_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_k_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )
    if ensemble_accepted_history.shape[1] != n_ensembles:
        raise ValueError(f"'ensemble_accepted_history' has {ensemble_accepted_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
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

    # outputs and work arrays, which the caller never sees
    super_ensembles = np.empty((max_group_size, n_ensembles*(n_ensembles-1),), dtype=np.int32, order='F')
    n_super_ensembles = ctypes.c_int(0)
    super_ensembles_overlap_coefficient = np.empty((max_group_size-1, n_ensembles*(n_ensembles-1),), dtype=np.float64, order='F')
    eligible = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_stop_condition = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_dimension = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    eligible_by_var_explained = np.empty((n_ensembles,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.ensemble_reconciliation_c(
        ensemble_masks,
        ensemble_stop_reason,
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        ensemble_U_history,
        ensemble_d_history,
        ensemble_S_history,
        ensemble_mu_history,
        ensemble_G_history,
        ensemble_k_history,
        ensemble_accepted_history,
        ctypes.byref(ctypes.c_int(o)),
        mode,
        ctypes.byref(ctypes.c_double(min_overlap_coefficient)),
        ctypes.byref(ctypes.c_bool(report_overlap_coefficient)),
        allowed_stop_reasons,
        None if filter_dim_min is None else ctypes.byref(ctypes.c_int(filter_dim_min)),
        None if filter_dim_max is None else ctypes.byref(ctypes.c_int(filter_dim_max)),
        None if var_explained_min is None else ctypes.byref(ctypes.c_double(var_explained_min)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        super_ensembles,
        ctypes.byref(n_super_ensembles),
        super_ensembles_overlap_coefficient,
        eligible,
        eligible_by_stop_condition,
        eligible_by_dimension,
        eligible_by_var_explained,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ENSEMBLE_RECONCILIATION_ARGUMENTS, _ENSEMBLE_RECONCILIATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    super_ensembles.flags.writeable = False
    super_ensembles_overlap_coefficient.flags.writeable = False
    eligible.flags.writeable = False
    eligible_by_stop_condition.flags.writeable = False
    eligible_by_dimension.flags.writeable = False
    eligible_by_var_explained.flags.writeable = False

    return {
        "super_ensembles": super_ensembles,
        "n_super_ensembles": n_super_ensembles.value,
        "super_ensembles_overlap_coefficient": super_ensembles_overlap_coefficient,
        "eligible": eligible,
        "eligible_by_stop_condition": eligible_by_stop_condition,
        "eligible_by_dimension": eligible_by_dimension,
        "eligible_by_var_explained": eligible_by_var_explained,
    }

def merge_to_super_ensembles(
        ensemble_masks,
        eligible,
        max_group_size,
        mode='report',
        min_overlap_coefficient=0.9,
        report_overlap_coefficient=False,
):
    r"""Group/report intersections among eligible ensembles into super-ensembles

    Parameters
    ----------
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_ensembles,), column-major (order='F')
        Per-ensemble membership, see Ensemble Identification's merged output
    eligible : np.ndarray[np.bool_] of shape (n_ensembles,)
        Per-ensemble eligibility to contribute a pair here at all -- see
        `tox_shape_truthful_clustering_filter_impl`'s own `filter_ensembles_impl`,
        this kernel's own sibling in `ensemble_reconciliation`'s two-call orchestration
    mode : str, one of 'report' | 'merge_overlap_coefficient' | 'merge_any', optional, default 'report'
        How intersections are processed

        The default value is `'report'`.
    min_overlap_coefficient : float, optional, default 0.9
        Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
        \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
        ``'merge_overlap_coefficient'``;
        ignored in every other mode
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    report_overlap_coefficient : bool, optional, default False
        Whether to compute and return `super_ensembles_overlap_coefficient` at all --
        see the note above on this being guarded, not unconditional
        The default value is `False`.
    max_group_size : int
        Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
        can hold; sizes its row dimension. `misc/mod_STC.md` suggests
        $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
        optional with an auto-applied default here, for the same reason as
        `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
        possibly-absent optional dummy, and a runtime-dependent value like
        $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
        default would need to be either.
        The minimum valid value is `2`.
        The maximum valid value is `n_ensembles`.

    Returns
    -------
    dict
        with keys:

        super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_ensembles*(n_ensembles-1),), column-major (order='F'), read-only
            One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            the group's actual size, and 0 in every row of an unused trailing column beyond
            `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            ``'report'``'s
            own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            intersects) -- deliberately not divided by 2, see `ensemble_reconciliation_impl`'s
            own identical note.
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_super_ensembles : int
            Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            actually filled
        super_ensembles_overlap_coefficient : np.ndarray[np.float64] of shape (max_group_size-1, n_ensembles*(n_ensembles-1),), column-major (order='F'), read-only
            Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            `report_overlap_coefficient` was requested -- see the note above.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_reconciliation::merge_to_super_ensembles`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ensemble_masks = np.asfortranarray(ensemble_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_masks' must be an array of np.bool_: {error}") from None
    if ensemble_masks.ndim != 2:
        raise ValueError(f"'ensemble_masks' must have 2 dimensions, but has {ensemble_masks.ndim}")
    try:
        eligible = np.ascontiguousarray(eligible, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eligible' must be an array of np.bool_: {error}") from None
    if eligible.ndim != 1:
        raise ValueError(f"'eligible' must have 1 dimension, but has {eligible.ndim}")
    mode = np.array([str(mode).lower().encode().ljust(25)], dtype="S25")

    # what the inputs already say, rather than asking for it again
    n_vectors = ensemble_masks.shape[0]
    n_ensembles = ensemble_masks.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if eligible.shape[0] != n_ensembles:
        raise ValueError(f"'eligible' has {eligible.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_ensembles == {n_ensembles}"
        )

    # outputs and work arrays, which the caller never sees
    super_ensembles = np.empty((max_group_size, n_ensembles*(n_ensembles-1),), dtype=np.int32, order='F')
    n_super_ensembles = ctypes.c_int(0)
    super_ensembles_overlap_coefficient = np.empty((max_group_size-1, n_ensembles*(n_ensembles-1),), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.merge_to_super_ensembles_c(
        ensemble_masks,
        eligible,
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        mode,
        ctypes.byref(ctypes.c_double(min_overlap_coefficient)),
        ctypes.byref(ctypes.c_bool(report_overlap_coefficient)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        super_ensembles,
        ctypes.byref(n_super_ensembles),
        super_ensembles_overlap_coefficient,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MERGE_TO_SUPER_ENSEMBLES_ARGUMENTS, _MERGE_TO_SUPER_ENSEMBLES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    super_ensembles.flags.writeable = False
    super_ensembles_overlap_coefficient.flags.writeable = False

    return {
        "super_ensembles": super_ensembles,
        "n_super_ensembles": n_super_ensembles.value,
        "super_ensembles_overlap_coefficient": super_ensembles_overlap_coefficient,
    }
