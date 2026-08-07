"""tox_shape_truthful_clustering_reconciliation

Generated from the kernel; do not edit -- regenerate instead.

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
_ENSEMBLE_RECONCILIATION_ARGUMENTS = ("ensemble_masks", "n_vectors", "n_ensembles", "mode", "min_jsi", "report_jsi", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_JSI", "ierr",)
#: For a derived argument, the one the caller passed it in
_ENSEMBLE_RECONCILIATION_ARGUMENT_SOURCES = (None, "ensemble_masks", "ensemble_masks", None, None, None, "super_ensembles", None, None, None, None,)

def ensemble_reconciliation(
        ensemble_masks,
        max_group_size,
        mode='report',
        min_jsi=0.1,
        report_jsi=False,
):
    r"""Identify and group intersecting ensembles from Ensemble Identification's merged output

    Parameters
    ----------
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_ensembles,), column-major (order='F')
        Per-ensemble membership, see Ensemble Identification's merged output
    mode : str, one of 'report' | 'merge_jsi' | 'merge_any', optional, default 'report'
        How intersections are processed

        The default value is `1`.
    min_jsi : float, optional, default 0.1
        Minimum Jaccard Similarity Index for an edge to qualify in mode
        ``'merge_jsi'``;
        ignored in every other mode
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    report_jsi : bool, optional, default False
        Whether to compute and return `super_ensembles_JSI` at all -- see the note
        above on this being guarded, not unconditional
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
            intersects) -- deliberately not divided by 2: the generator translates this
            specification expression close to verbatim into the Python/R bindings, where
            `/` on two integers is true division, not Fortran's own truncating integer
            division, so a literal `/2` here breaks the generated Python binding (a `float`
            where `np.empty`'s shape wants an `int`); see `misc/code_gen_footgun.md`. A
            safe, if looser, upper bound for modes 2 and 3 too, whose groups can never
            outnumber mode 1's own worst case.
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_super_ensembles : int
            Number of leading columns of `super_ensembles`/`super_ensembles_JSI` actually
            filled
        super_ensembles_JSI : np.ndarray[np.float64] of shape (max_group_size-1, n_ensembles*(n_ensembles-1),), column-major (order='F'), read-only
            Column $l$, row $c_i$: the JSI between the ensembles in `super_ensembles(c_i, l)`
            and `super_ensembles(c_i + 1, l)`. All zero unless `report_jsi` was requested --
            see the note above.
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
    mode = np.array([str(mode).lower().encode()], dtype="S9")

    # what the inputs already say, rather than asking for it again
    n_vectors = ensemble_masks.shape[0]
    n_ensembles = ensemble_masks.shape[1]

    # outputs and work arrays, which the caller never sees
    super_ensembles = np.empty((max_group_size, n_ensembles*(n_ensembles-1),), dtype=np.int32, order='F')
    n_super_ensembles = ctypes.c_int(0)
    super_ensembles_JSI = np.empty((max_group_size-1, n_ensembles*(n_ensembles-1),), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.ensemble_reconciliation_c(
        ensemble_masks,
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_ensembles)),
        mode,
        ctypes.byref(ctypes.c_double(min_jsi)),
        ctypes.byref(ctypes.c_bool(report_jsi)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        super_ensembles,
        ctypes.byref(n_super_ensembles),
        super_ensembles_JSI,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ENSEMBLE_RECONCILIATION_ARGUMENTS, _ENSEMBLE_RECONCILIATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    super_ensembles.flags.writeable = False
    super_ensembles_JSI.flags.writeable = False

    return {
        "super_ensembles": super_ensembles,
        "n_super_ensembles": n_super_ensembles.value,
        "super_ensembles_JSI": super_ensembles_JSI,
    }
