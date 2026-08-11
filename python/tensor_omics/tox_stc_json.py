"""tox_stc_json

Serialization of Shape Truthful Clustering (STC) pipeline results into the JSON format

Python binding, generated from tox_stc_json. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.serialize_stc_results_as_json_c.restype = None
_lib.serialize_stc_results_as_json_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
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
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_STC_RESULTS_AS_JSON_ARGUMENTS = ("filename", "n_dimensions", "n_vectors", "n_selected_seed", "o", "max_group_size", "n_super_ensembles", "vectors", "dim_names", "seed_selection_mask", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_low_confidence_masks", "super_ensembles", "k_min", "k_density", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "exclusion_radius_percentile", "bandwidth_percentile", "reconciliation_mode", "min_overlap_coefficient", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_STC_RESULTS_AS_JSON_ARGUMENT_SOURCES = (None, "vectors", "vectors", "ensemble_masks", "ensemble_U_history", "super_ensembles", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.write_stc_interactive_html_report_c.restype = None
_lib.write_stc_interactive_html_report_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
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
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_double)),
    nullable(ctypes.POINTER(ctypes.c_int)),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_WRITE_STC_INTERACTIVE_HTML_REPORT_ARGUMENTS = ("filename", "n_dimensions", "n_vectors", "n_selected_seed", "o", "max_group_size", "n_super_ensembles", "vectors", "dim_names", "seed_selection_mask", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_low_confidence_masks", "super_ensembles", "k_min", "k_density", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "f_max", "a", "exclusion_radius_percentile", "bandwidth_percentile", "reconciliation_mode", "min_overlap_coefficient", "estimated_k_min", "estimated_k_density", "estimated_density_quantile", "estimated_chordal_dist_max_as_prcnt_of_range", "estimated_G_max", "estimated_d_max", "ierr",)
#: For a derived argument, the one the caller passed it in
_WRITE_STC_INTERACTIVE_HTML_REPORT_ARGUMENT_SOURCES = (None, "vectors", "vectors", "ensemble_masks", "ensemble_U_history", "super_ensembles", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def serialize_stc_results_as_json(
        filename,
        n_super_ensembles,
        vectors,
        dim_names,
        seed_selection_mask,
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_growth_radii,
        ensemble_U_history,
        ensemble_S_history,
        ensemble_d_history,
        ensemble_G_history,
        ensemble_mu_history,
        ensemble_k_history,
        ensemble_low_confidence_masks,
        super_ensembles,
        k_min,
        k_density,
        chordal_dist_max_as_prcnt_of_range,
        d_max,
        G_max,
        RMSE_change_max,
        f_max,
        a,
        exclusion_radius_percentile,
        bandwidth_percentile,
        reconciliation_mode,
        min_overlap_coefficient,
        estimated_k_min=None,
        estimated_k_density=None,
        estimated_density_quantile=None,
        estimated_chordal_dist_max_as_prcnt_of_range=None,
        estimated_G_max=None,
        estimated_d_max=None,
):
    r"""Serializes an STC run's raw pipeline results as JSON

    Parameters
    ----------
    filename : str
        Name of the JSON file to write
    n_super_ensembles : int
        Number of leading columns of `super_ensembles` actually filled
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    dim_names : sequence of str, of length n_dimensions
        Per-dimension display name
    seed_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Seed selection, see `seeds`
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble accepted membership, one column per seed
    ensemble_stop_reason : np.ndarray[np.int32] of shape (n_selected_seed,)
        Per-ensemble Stop Condition
    ensemble_growth_radii : np.ndarray[np.float64] of shape (n_selected_seed,)
        Per-ensemble growth radius
    ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing tangent+normal bases
    ensemble_S_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing singular values
    ensemble_d_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing intrinsic dimensions
    ensemble_G_history : np.ndarray[np.float64] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing spectral gaps
    ensemble_mu_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing centers
    ensemble_k_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing sizes
    ensemble_low_confidence_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble iteration-1 fallback membership
    super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_selected_seed*(n_selected_seed-1),), column-major (order='F')
        One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
    k_min : int
        This run's neighborhood size for each seed's growth radius
    k_density : int
        This run's density estimation neighborhood size
    chordal_dist_max_as_prcnt_of_range : float
        This run's maximum tolerated chordal distance between tangent bases
    d_max : int
        This run's maximum tolerated change in intrinsic dimension
    G_max : float
        This run's maximum tolerated |log(G_tp1/G_t)|
    RMSE_change_max : float
        This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
    f_max : float
        This run's ensemble size fraction of N above which growth is abandoned
    a : int
        This run's minimum accepted-iteration count for a stable rejection
    exclusion_radius_percentile : float
        This run's seeding exclusion radius percentile
    bandwidth_percentile : float
        This run's density-estimate kernel bandwidth percentile
    reconciliation_mode : str, one of 'report' | 'merge_overlap_coefficient' | 'merge_any'
        This run's `ensemble_reconciliation` mode

    min_overlap_coefficient : float
        This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
    estimated_k_min : int, optional
        `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
    estimated_k_density : int, optional
        `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
    estimated_density_quantile : float, optional
        `estimate_stc_parameters`'s proposed density quantile, if estimation was used
    estimated_chordal_dist_max_as_prcnt_of_range : float, optional
        `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
        estimation was used
    estimated_G_max : float, optional
        `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
    estimated_d_max : int, optional
        `estimate_stc_parameters`'s proposed `d_max`, if estimation was used

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc_json::serialize_stc_results_as_json`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        vectors = np.asfortranarray(vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors' must be an array of np.float64: {error}") from None
    if vectors.ndim != 2:
        raise ValueError(f"'vectors' must have 2 dimensions, but has {vectors.ndim}")
    try:
        dim_names = np.asarray([str(_s).encode() for _s in dim_names], dtype="S")
    except TypeError as error:
        raise TypeError(f"'dim_names' must be a sequence of strings: {error}") from None
    if dim_names.ndim != 1:
        raise ValueError(f"'dim_names' must have 1 dimension, but has {dim_names.ndim}")
    try:
        seed_selection_mask = np.ascontiguousarray(seed_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'seed_selection_mask' must be an array of np.bool_: {error}") from None
    if seed_selection_mask.ndim != 1:
        raise ValueError(f"'seed_selection_mask' must have 1 dimension, but has {seed_selection_mask.ndim}")
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
        ensemble_growth_radii = np.ascontiguousarray(ensemble_growth_radii, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_growth_radii' must be an array of np.float64: {error}") from None
    if ensemble_growth_radii.ndim != 1:
        raise ValueError(f"'ensemble_growth_radii' must have 1 dimension, but has {ensemble_growth_radii.ndim}")
    try:
        ensemble_U_history = np.asfortranarray(ensemble_U_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_U_history' must be an array of np.float64: {error}") from None
    if ensemble_U_history.ndim != 4:
        raise ValueError(f"'ensemble_U_history' must have 4 dimensions, but has {ensemble_U_history.ndim}")
    try:
        ensemble_S_history = np.asfortranarray(ensemble_S_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_S_history' must be an array of np.float64: {error}") from None
    if ensemble_S_history.ndim != 3:
        raise ValueError(f"'ensemble_S_history' must have 3 dimensions, but has {ensemble_S_history.ndim}")
    try:
        ensemble_d_history = np.asfortranarray(ensemble_d_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_d_history' must be an array of np.int32: {error}") from None
    if ensemble_d_history.ndim != 2:
        raise ValueError(f"'ensemble_d_history' must have 2 dimensions, but has {ensemble_d_history.ndim}")
    try:
        ensemble_G_history = np.asfortranarray(ensemble_G_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_G_history' must be an array of np.float64: {error}") from None
    if ensemble_G_history.ndim != 2:
        raise ValueError(f"'ensemble_G_history' must have 2 dimensions, but has {ensemble_G_history.ndim}")
    try:
        ensemble_mu_history = np.asfortranarray(ensemble_mu_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_mu_history' must be an array of np.float64: {error}") from None
    if ensemble_mu_history.ndim != 3:
        raise ValueError(f"'ensemble_mu_history' must have 3 dimensions, but has {ensemble_mu_history.ndim}")
    try:
        ensemble_k_history = np.asfortranarray(ensemble_k_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_k_history' must be an array of np.int32: {error}") from None
    if ensemble_k_history.ndim != 2:
        raise ValueError(f"'ensemble_k_history' must have 2 dimensions, but has {ensemble_k_history.ndim}")
    try:
        ensemble_low_confidence_masks = np.asfortranarray(ensemble_low_confidence_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_low_confidence_masks' must be an array of np.bool_: {error}") from None
    if ensemble_low_confidence_masks.ndim != 2:
        raise ValueError(f"'ensemble_low_confidence_masks' must have 2 dimensions, but has {ensemble_low_confidence_masks.ndim}")
    try:
        super_ensembles = np.asfortranarray(super_ensembles, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'super_ensembles' must be an array of np.int32: {error}") from None
    if super_ensembles.ndim != 2:
        raise ValueError(f"'super_ensembles' must have 2 dimensions, but has {super_ensembles.ndim}")
    reconciliation_mode = np.array([str(reconciliation_mode).lower().encode()], dtype="S25")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_selected_seed = int(seed_selection_mask.sum())
    o = ensemble_U_history.shape[2]
    max_group_size = super_ensembles.shape[0]
    dim_names_strlen = dim_names.itemsize

    # Fortran cannot check that shared extents agree; this can
    if dim_names.shape[0] != n_dimensions:
        raise ValueError(f"'dim_names' has {dim_names.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_U_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_U_history' has {ensemble_U_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_S_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_mu_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if seed_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'seed_selection_mask' has {seed_selection_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_masks' has {ensemble_masks.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_low_confidence_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_stop_reason.shape[0] != n_selected_seed:
        raise ValueError(f"'ensemble_stop_reason' has {ensemble_stop_reason.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_growth_radii.shape[0] != n_selected_seed:
        raise ValueError(f"'ensemble_growth_radii' has {ensemble_growth_radii.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_U_history.shape[3] != n_selected_seed:
        raise ValueError(f"'ensemble_U_history' has {ensemble_U_history.shape[3]} along axis 3, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_S_history.shape[2] != n_selected_seed:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_d_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_G_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_mu_history.shape[2] != n_selected_seed:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_k_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_low_confidence_masks.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_S_history.shape[1] != o:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_d_history.shape[0] != o:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_G_history.shape[0] != o:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_mu_history.shape[1] != o:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_k_history.shape[0] != o:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_stc_results_as_json_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_selected_seed)),
        ctypes.byref(ctypes.c_int(o)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        ctypes.byref(ctypes.c_int(n_super_ensembles)),
        vectors,
        dim_names,
        ctypes.byref(ctypes.c_int(dim_names_strlen)),
        seed_selection_mask,
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_growth_radii,
        ensemble_U_history,
        ensemble_S_history,
        ensemble_d_history,
        ensemble_G_history,
        ensemble_mu_history,
        ensemble_k_history,
        ensemble_low_confidence_masks,
        super_ensembles,
        ctypes.byref(ctypes.c_int(k_min)),
        ctypes.byref(ctypes.c_int(k_density)),
        ctypes.byref(ctypes.c_double(chordal_dist_max_as_prcnt_of_range)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(ctypes.c_double(RMSE_change_max)),
        ctypes.byref(ctypes.c_double(f_max)),
        ctypes.byref(ctypes.c_int(a)),
        ctypes.byref(ctypes.c_double(exclusion_radius_percentile)),
        ctypes.byref(ctypes.c_double(bandwidth_percentile)),
        reconciliation_mode,
        ctypes.byref(ctypes.c_double(min_overlap_coefficient)),
        None if estimated_k_min is None else ctypes.byref(ctypes.c_int(estimated_k_min)),
        None if estimated_k_density is None else ctypes.byref(ctypes.c_int(estimated_k_density)),
        None if estimated_density_quantile is None else ctypes.byref(ctypes.c_double(estimated_density_quantile)),
        None if estimated_chordal_dist_max_as_prcnt_of_range is None else ctypes.byref(ctypes.c_double(estimated_chordal_dist_max_as_prcnt_of_range)),
        None if estimated_G_max is None else ctypes.byref(ctypes.c_double(estimated_G_max)),
        None if estimated_d_max is None else ctypes.byref(ctypes.c_int(estimated_d_max)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_STC_RESULTS_AS_JSON_ARGUMENTS, _SERIALIZE_STC_RESULTS_AS_JSON_ARGUMENT_SOURCES)

    return None

def write_stc_interactive_html_report(
        filename,
        n_super_ensembles,
        vectors,
        dim_names,
        seed_selection_mask,
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_growth_radii,
        ensemble_U_history,
        ensemble_S_history,
        ensemble_d_history,
        ensemble_G_history,
        ensemble_mu_history,
        ensemble_k_history,
        ensemble_low_confidence_masks,
        super_ensembles,
        k_min,
        k_density,
        chordal_dist_max_as_prcnt_of_range,
        d_max,
        G_max,
        RMSE_change_max,
        f_max,
        a,
        exclusion_radius_percentile,
        bandwidth_percentile,
        reconciliation_mode,
        min_overlap_coefficient,
        estimated_k_min=None,
        estimated_k_density=None,
        estimated_density_quantile=None,
        estimated_chordal_dist_max_as_prcnt_of_range=None,
        estimated_G_max=None,
        estimated_d_max=None,
):
    r"""Writes a self-contained interactive HTML report for an STC run

    Parameters
    ----------
    filename : str
        Name of the HTML file to write
    n_super_ensembles : int
        Number of leading columns of `super_ensembles` actually filled
    vectors : np.ndarray[np.float64] of shape (n_dimensions, n_vectors,), column-major (order='F')
        Input data matrix
    dim_names : sequence of str, of length n_dimensions
        Per-dimension display name
    seed_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Seed selection, see `seeds`
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble accepted membership, one column per seed
    ensemble_stop_reason : np.ndarray[np.int32] of shape (n_selected_seed,)
        Per-ensemble Stop Condition
    ensemble_growth_radii : np.ndarray[np.float64] of shape (n_selected_seed,)
        Per-ensemble growth radius
    ensemble_U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing tangent+normal bases
    ensemble_S_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing singular values
    ensemble_d_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing intrinsic dimensions
    ensemble_G_history : np.ndarray[np.float64] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing spectral gaps
    ensemble_mu_history : np.ndarray[np.float64] of shape (n_dimensions, o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing centers
    ensemble_k_history : np.ndarray[np.int32] of shape (o, n_selected_seed,), column-major (order='F')
        Per-ensemble trailing sizes
    ensemble_low_confidence_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble iteration-1 fallback membership
    super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_selected_seed*(n_selected_seed-1),), column-major (order='F')
        One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
    k_min : int
        This run's neighborhood size for each seed's growth radius
    k_density : int
        This run's density estimation neighborhood size
    chordal_dist_max_as_prcnt_of_range : float
        This run's maximum tolerated chordal distance between tangent bases
    d_max : int
        This run's maximum tolerated change in intrinsic dimension
    G_max : float
        This run's maximum tolerated |log(G_tp1/G_t)|
    RMSE_change_max : float
        This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
    f_max : float
        This run's ensemble size fraction of N above which growth is abandoned
    a : int
        This run's minimum accepted-iteration count for a stable rejection
    exclusion_radius_percentile : float
        This run's seeding exclusion radius percentile
    bandwidth_percentile : float
        This run's density-estimate kernel bandwidth percentile
    reconciliation_mode : str, one of 'report' | 'merge_overlap_coefficient' | 'merge_any'
        This run's `ensemble_reconciliation` mode

    min_overlap_coefficient : float
        This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
    estimated_k_min : int, optional
        `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
    estimated_k_density : int, optional
        `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
    estimated_density_quantile : float, optional
        `estimate_stc_parameters`'s proposed density quantile, if estimation was used
    estimated_chordal_dist_max_as_prcnt_of_range : float, optional
        `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
        estimation was used
    estimated_G_max : float, optional
        `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
    estimated_d_max : int, optional
        `estimate_stc_parameters`'s proposed `d_max`, if estimation was used

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc_json::write_stc_interactive_html_report`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        vectors = np.asfortranarray(vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors' must be an array of np.float64: {error}") from None
    if vectors.ndim != 2:
        raise ValueError(f"'vectors' must have 2 dimensions, but has {vectors.ndim}")
    try:
        dim_names = np.asarray([str(_s).encode() for _s in dim_names], dtype="S")
    except TypeError as error:
        raise TypeError(f"'dim_names' must be a sequence of strings: {error}") from None
    if dim_names.ndim != 1:
        raise ValueError(f"'dim_names' must have 1 dimension, but has {dim_names.ndim}")
    try:
        seed_selection_mask = np.ascontiguousarray(seed_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'seed_selection_mask' must be an array of np.bool_: {error}") from None
    if seed_selection_mask.ndim != 1:
        raise ValueError(f"'seed_selection_mask' must have 1 dimension, but has {seed_selection_mask.ndim}")
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
        ensemble_growth_radii = np.ascontiguousarray(ensemble_growth_radii, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_growth_radii' must be an array of np.float64: {error}") from None
    if ensemble_growth_radii.ndim != 1:
        raise ValueError(f"'ensemble_growth_radii' must have 1 dimension, but has {ensemble_growth_radii.ndim}")
    try:
        ensemble_U_history = np.asfortranarray(ensemble_U_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_U_history' must be an array of np.float64: {error}") from None
    if ensemble_U_history.ndim != 4:
        raise ValueError(f"'ensemble_U_history' must have 4 dimensions, but has {ensemble_U_history.ndim}")
    try:
        ensemble_S_history = np.asfortranarray(ensemble_S_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_S_history' must be an array of np.float64: {error}") from None
    if ensemble_S_history.ndim != 3:
        raise ValueError(f"'ensemble_S_history' must have 3 dimensions, but has {ensemble_S_history.ndim}")
    try:
        ensemble_d_history = np.asfortranarray(ensemble_d_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_d_history' must be an array of np.int32: {error}") from None
    if ensemble_d_history.ndim != 2:
        raise ValueError(f"'ensemble_d_history' must have 2 dimensions, but has {ensemble_d_history.ndim}")
    try:
        ensemble_G_history = np.asfortranarray(ensemble_G_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_G_history' must be an array of np.float64: {error}") from None
    if ensemble_G_history.ndim != 2:
        raise ValueError(f"'ensemble_G_history' must have 2 dimensions, but has {ensemble_G_history.ndim}")
    try:
        ensemble_mu_history = np.asfortranarray(ensemble_mu_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_mu_history' must be an array of np.float64: {error}") from None
    if ensemble_mu_history.ndim != 3:
        raise ValueError(f"'ensemble_mu_history' must have 3 dimensions, but has {ensemble_mu_history.ndim}")
    try:
        ensemble_k_history = np.asfortranarray(ensemble_k_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_k_history' must be an array of np.int32: {error}") from None
    if ensemble_k_history.ndim != 2:
        raise ValueError(f"'ensemble_k_history' must have 2 dimensions, but has {ensemble_k_history.ndim}")
    try:
        ensemble_low_confidence_masks = np.asfortranarray(ensemble_low_confidence_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_low_confidence_masks' must be an array of np.bool_: {error}") from None
    if ensemble_low_confidence_masks.ndim != 2:
        raise ValueError(f"'ensemble_low_confidence_masks' must have 2 dimensions, but has {ensemble_low_confidence_masks.ndim}")
    try:
        super_ensembles = np.asfortranarray(super_ensembles, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'super_ensembles' must be an array of np.int32: {error}") from None
    if super_ensembles.ndim != 2:
        raise ValueError(f"'super_ensembles' must have 2 dimensions, but has {super_ensembles.ndim}")
    reconciliation_mode = np.array([str(reconciliation_mode).lower().encode()], dtype="S25")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    n_dimensions = vectors.shape[0]
    n_vectors = vectors.shape[1]
    n_selected_seed = int(seed_selection_mask.sum())
    o = ensemble_U_history.shape[2]
    max_group_size = super_ensembles.shape[0]
    dim_names_strlen = dim_names.itemsize

    # Fortran cannot check that shared extents agree; this can
    if dim_names.shape[0] != n_dimensions:
        raise ValueError(f"'dim_names' has {dim_names.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_U_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_U_history' has {ensemble_U_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_S_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if ensemble_mu_history.shape[0] != n_dimensions:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[0]} along axis 0, but "
            f"'vectors' implies n_dimensions == {n_dimensions}"
        )
    if seed_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'seed_selection_mask' has {seed_selection_mask.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_masks' has {ensemble_masks.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_low_confidence_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[0]} along axis 0, but "
            f"'vectors' implies n_vectors == {n_vectors}"
        )
    if ensemble_stop_reason.shape[0] != n_selected_seed:
        raise ValueError(f"'ensemble_stop_reason' has {ensemble_stop_reason.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_growth_radii.shape[0] != n_selected_seed:
        raise ValueError(f"'ensemble_growth_radii' has {ensemble_growth_radii.shape[0]} along axis 0, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_U_history.shape[3] != n_selected_seed:
        raise ValueError(f"'ensemble_U_history' has {ensemble_U_history.shape[3]} along axis 3, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_S_history.shape[2] != n_selected_seed:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_d_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_G_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_mu_history.shape[2] != n_selected_seed:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[2]} along axis 2, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_k_history.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_low_confidence_masks.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )
    if ensemble_S_history.shape[1] != o:
        raise ValueError(f"'ensemble_S_history' has {ensemble_S_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_d_history.shape[0] != o:
        raise ValueError(f"'ensemble_d_history' has {ensemble_d_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_G_history.shape[0] != o:
        raise ValueError(f"'ensemble_G_history' has {ensemble_G_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_mu_history.shape[1] != o:
        raise ValueError(f"'ensemble_mu_history' has {ensemble_mu_history.shape[1]} along axis 1, but "
            f"'ensemble_U_history' implies o == {o}"
        )
    if ensemble_k_history.shape[0] != o:
        raise ValueError(f"'ensemble_k_history' has {ensemble_k_history.shape[0]} along axis 0, but "
            f"'ensemble_U_history' implies o == {o}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.write_stc_interactive_html_report_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_selected_seed)),
        ctypes.byref(ctypes.c_int(o)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        ctypes.byref(ctypes.c_int(n_super_ensembles)),
        vectors,
        dim_names,
        ctypes.byref(ctypes.c_int(dim_names_strlen)),
        seed_selection_mask,
        ensemble_masks,
        ensemble_stop_reason,
        ensemble_growth_radii,
        ensemble_U_history,
        ensemble_S_history,
        ensemble_d_history,
        ensemble_G_history,
        ensemble_mu_history,
        ensemble_k_history,
        ensemble_low_confidence_masks,
        super_ensembles,
        ctypes.byref(ctypes.c_int(k_min)),
        ctypes.byref(ctypes.c_int(k_density)),
        ctypes.byref(ctypes.c_double(chordal_dist_max_as_prcnt_of_range)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(ctypes.c_double(RMSE_change_max)),
        ctypes.byref(ctypes.c_double(f_max)),
        ctypes.byref(ctypes.c_int(a)),
        ctypes.byref(ctypes.c_double(exclusion_radius_percentile)),
        ctypes.byref(ctypes.c_double(bandwidth_percentile)),
        reconciliation_mode,
        ctypes.byref(ctypes.c_double(min_overlap_coefficient)),
        None if estimated_k_min is None else ctypes.byref(ctypes.c_int(estimated_k_min)),
        None if estimated_k_density is None else ctypes.byref(ctypes.c_int(estimated_k_density)),
        None if estimated_density_quantile is None else ctypes.byref(ctypes.c_double(estimated_density_quantile)),
        None if estimated_chordal_dist_max_as_prcnt_of_range is None else ctypes.byref(ctypes.c_double(estimated_chordal_dist_max_as_prcnt_of_range)),
        None if estimated_G_max is None else ctypes.byref(ctypes.c_double(estimated_G_max)),
        None if estimated_d_max is None else ctypes.byref(ctypes.c_int(estimated_d_max)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _WRITE_STC_INTERACTIVE_HTML_REPORT_ARGUMENTS, _WRITE_STC_INTERACTIVE_HTML_REPORT_ARGUMENT_SOURCES)

    return None
