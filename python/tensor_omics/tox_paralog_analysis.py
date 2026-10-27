r"""tox_paralog_analysis

Paralog-subset expression patterns -- dosage effect and subfunctionalization -- relative to an ancestral ortholog.

Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
starting from single genes. At every extension step the candidate is scored against the pattern-specific
criterion (small angle plus magnitude gain for dosage effect, or bounded residual distance for
subfunctionalization); subsets that can no longer satisfy the criterion are pruned instead of being
extended further, which keeps the combinatorial subset search tractable.

There is a routine per pattern rather than one taking a pattern flag:
`detect_dosage_effect` and `detect_subfunctionalization` find the subsets, and
`filter_paralogs_by_pattern_dosage_effect` / `_subfunctionalization` reduce a set of genes to
those matching. Which criterion a result was scored against is therefore visible at the call
site, and each routine takes only the thresholds its own pattern needs.

Python binding, generated from tox_paralog_analysis. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.detect_neofunctionalization_c.restype = None
_lib.detect_neofunctionalization_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_NEOFUNCTIONALIZATION_ARGUMENTS = ("ancestors", "n_families", "genes", "n_axes", "gene_to_fam", "n_genes", "thresholds", "neofunc", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETECT_NEOFUNCTIONALIZATION_ARGUMENT_SOURCES = (None, "ancestors", None, "ancestors", None, "genes", None, None, None,)

_lib.detect_dosage_effect_c.restype = None
_lib.detect_dosage_effect_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_DOSAGE_EFFECT_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "max_angle", "gain_gamma", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETECT_DOSAGE_EFFECT_ARGUMENT_SOURCES = (None, None, "genes", "ancestor", None, "filtered_paralogs_mask", None, None, None, "work_arr_paralog_subsets", None, None, None,)

_lib.detect_subfunctionalization_c.restype = None
_lib.detect_subfunctionalization_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_SUBFUNCTIONALIZATION_ARGUMENTS = ("ancestor", "genes", "n_genes", "n_dims", "filtered_paralogs_mask", "n_mask_chunks", "n_results", "max_subset_size", "work_arr_paralog_subsets", "n_paralog_subsets", "rdi_threshold", "paralog_norms", "sorted_paralog_norms_perm", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETECT_SUBFUNCTIONALIZATION_ARGUMENT_SOURCES = (None, None, "genes", "ancestor", None, "filtered_paralogs_mask", None, None, None, "work_arr_paralog_subsets", None, None, None, None,)

_lib.filter_paralogs_by_pattern_dosage_effect_c.restype = None
_lib.filter_paralogs_by_pattern_dosage_effect_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENTS = ("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENT_SOURCES = (None, None, "gene_angles", "masks", None, None, "masks", None,)

_lib.filter_paralogs_by_pattern_subfunctionalization_c.restype = None
_lib.filter_paralogs_by_pattern_subfunctionalization_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENTS = ("gene_angles", "threshold", "n_genes", "n_families", "gene_to_fam", "masks", "n_mask_chunks", "ierr",)
#: For a derived argument, the one the caller passed it in
_FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENT_SOURCES = (None, None, "gene_angles", "masks", None, None, "masks", None,)

def detect_neofunctionalization(
        ancestors,
        genes,
        gene_to_fam,
        thresholds,
):
    r"""Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis

    Parameters
    ----------
    ancestors : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        RAP projected unit length expression vector of ancestral ortholog
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.
    genes : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        RAP projected unit length expression vectors of genes
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    thresholds : np.ndarray[np.float64] of shape (n_axes,)
        threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    neofunc : np.ndarray[np.bool_] of shape (n_genes, n_axes,), column-major (order='F'), read-only
        `True` if neofunctionalization has been detected for the respective axes, always `False` for unassigned genes
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
    Generated from the Fortran procedure `tox_paralog_analysis::detect_neofunctionalization`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestors = _np.asfortranarray(ancestors, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'ancestors' must be an array of np.float64: {_error}", argument="ancestors") from None
    if ancestors.ndim != 2:
        raise _ToxArgumentError(f"'ancestors' must have 2 dimensions, but has {ancestors.ndim}", argument="ancestors")
    try:
        genes = _np.asfortranarray(genes, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'genes' must be an array of np.float64: {_error}", argument="genes") from None
    if genes.ndim != 2:
        raise _ToxArgumentError(f"'genes' must have 2 dimensions, but has {genes.ndim}", argument="genes")
    try:
        gene_to_fam = _np.ascontiguousarray(gene_to_fam, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_fam' must be an array of np.int32: {_error}", argument="gene_to_fam") from None
    if gene_to_fam.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}", argument="gene_to_fam")
    try:
        thresholds = _np.ascontiguousarray(thresholds, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'thresholds' must be an array of np.float64: {_error}", argument="thresholds") from None
    if thresholds.ndim != 1:
        raise _ToxArgumentError(f"'thresholds' must have 1 dimension, but has {thresholds.ndim}", argument="thresholds")

    # what the inputs already say, rather than asking for it again
    n_families = ancestors.shape[1]
    n_axes = ancestors.shape[0]
    n_genes = genes.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if genes.shape[0] != n_axes:
        raise _ToxArgumentError(f"'genes' has {genes.shape[0]} along axis 0, but 'ancestors' implies n_axes == {n_axes}", argument="genes")
    if thresholds.shape[0] != n_axes:
        raise _ToxArgumentError(f"'thresholds' has {thresholds.shape[0]} along axis 0, but 'ancestors' implies n_axes == {n_axes}", argument="thresholds")
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'genes' implies n_genes == {n_genes}", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    neofunc = _np.empty((n_genes, n_axes,), dtype=_np.bool_, order='F')
    ierr = _ctypes.c_int(0)

    _lib.detect_neofunctionalization_c(
        ancestors,
        _ctypes.byref(_ctypes.c_int(n_families)),
        genes,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        gene_to_fam,
        _ctypes.byref(_ctypes.c_int(n_genes)),
        thresholds,
        neofunc,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETECT_NEOFUNCTIONALIZATION_ARGUMENTS, _DETECT_NEOFUNCTIONALIZATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neofunc.flags.writeable = False

    return neofunc

def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        max_angle=3.141592653589793,
        gain_gamma=0.1,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        :func:`tensor_omics.calc_work_arr_paralog_subsets_size` first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    max_angle : float, optional, default 3.141592653589793
        maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
        The default value is `4.0*atan(1.0)`.
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    gain_gamma : float, optional, default 0.1
        positive magnitude gain for dosage effect
        The default value is `0.1`.
        The minimum valid value is `above(0.0)`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F'), read-only
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks
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
    Generated from the Fortran procedure `tox_paralog_analysis::detect_dosage_effect`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestor = _np.ascontiguousarray(ancestor, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'ancestor' must be an array of np.float64: {_error}", argument="ancestor") from None
    if ancestor.ndim != 1:
        raise _ToxArgumentError(f"'ancestor' must have 1 dimension, but has {ancestor.ndim}", argument="ancestor")
    try:
        genes = _np.asfortranarray(genes, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'genes' must be an array of np.float64: {_error}", argument="genes") from None
    if genes.ndim != 2:
        raise _ToxArgumentError(f"'genes' must have 2 dimensions, but has {genes.ndim}", argument="genes")
    try:
        filtered_paralogs_mask = _np.ascontiguousarray(filtered_paralogs_mask, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'filtered_paralogs_mask' must be an array of np.int32: {_error}", argument="filtered_paralogs_mask") from None
    if filtered_paralogs_mask.ndim != 1:
        raise _ToxArgumentError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}", argument="filtered_paralogs_mask")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_impl import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

    # Fortran cannot check that shared extents agree; this can
    if genes.shape[0] != n_dims:
        raise _ToxArgumentError(f"'genes' has {genes.shape[0]} along axis 0, but 'ancestor' implies n_dims == {n_dims}", argument="genes")

    # outputs and work arrays, which the caller never sees
    n_results = _ctypes.c_int(0)
    work_arr_paralog_subsets = _np.empty((n_mask_chunks, n_paralog_subsets,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.detect_dosage_effect_c(
        ancestor,
        genes,
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        _ctypes.byref(_ctypes.c_int(n_mask_chunks)),
        _ctypes.byref(n_results),
        _ctypes.byref(_ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        _ctypes.byref(_ctypes.c_int(n_paralog_subsets)),
        _ctypes.byref(_ctypes.c_double(max_angle)),
        _ctypes.byref(_ctypes.c_double(gain_gamma)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETECT_DOSAGE_EFFECT_ARGUMENTS, _DETECT_DOSAGE_EFFECT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    work_arr_paralog_subsets.flags.writeable = False

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def detect_subfunctionalization(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        rdi_threshold,
        paralog_norms,
        sorted_paralog_norms_perm,
):
    r"""Identifies subsets of paralogs matching this pattern

    Parameters
    ----------
    ancestor : np.ndarray[np.float64] of shape (n_dims,)
        expression vector of ancestral ortholog
    genes : np.ndarray[np.float64] of shape (n_dims, n_genes,), column-major (order='F')
        expression vectors of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
    max_subset_size : int
        maximum subset size of checked gene subsets. Too large a value is capped to the
        maximum valid size. The bindings cap it automatically while sizing the work
        array; a Fortran caller caps it by calling
        :func:`tensor_omics.calc_work_arr_paralog_subsets_size` first.
        Zero is in range and means there is no subset to check -- the sizing routine reports
        it whenever the filtered families hold a single gene each. It reports a work array
        of zero slots along with it, which this routine does not accept, so a caller that
        gets zero back has nothing to detect and should not call here at all.
        The minimum valid value is `0`.
    rdi_threshold : float
        max allowed residual distance from `ancestor`
        The minimum valid value is `0.0`.
    paralog_norms : np.ndarray[np.float64] of shape (n_genes,)
        euclidean norms of the genes, used for subset pruning (`norm` from `f42_vector_impl` computes them)
        The minimum valid value is `0.0`.
    sorted_paralog_norms_perm : np.ndarray[np.int32] of shape (n_genes,)
        ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.

    Returns
    -------
    dict
        with keys:

        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        work_arr_paralog_subsets : np.ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets,), column-major (order='F'), read-only
            working array to hold bitmask encoded subsets for detection.
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks
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
    Generated from the Fortran procedure `tox_paralog_analysis::detect_subfunctionalization`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        ancestor = _np.ascontiguousarray(ancestor, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'ancestor' must be an array of np.float64: {_error}", argument="ancestor") from None
    if ancestor.ndim != 1:
        raise _ToxArgumentError(f"'ancestor' must have 1 dimension, but has {ancestor.ndim}", argument="ancestor")
    try:
        genes = _np.asfortranarray(genes, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'genes' must be an array of np.float64: {_error}", argument="genes") from None
    if genes.ndim != 2:
        raise _ToxArgumentError(f"'genes' must have 2 dimensions, but has {genes.ndim}", argument="genes")
    try:
        filtered_paralogs_mask = _np.ascontiguousarray(filtered_paralogs_mask, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'filtered_paralogs_mask' must be an array of np.int32: {_error}", argument="filtered_paralogs_mask") from None
    if filtered_paralogs_mask.ndim != 1:
        raise _ToxArgumentError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}", argument="filtered_paralogs_mask")
    try:
        paralog_norms = _np.ascontiguousarray(paralog_norms, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'paralog_norms' must be an array of np.float64: {_error}", argument="paralog_norms") from None
    if paralog_norms.ndim != 1:
        raise _ToxArgumentError(f"'paralog_norms' must have 1 dimension, but has {paralog_norms.ndim}", argument="paralog_norms")
    try:
        sorted_paralog_norms_perm = _np.ascontiguousarray(sorted_paralog_norms_perm, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'sorted_paralog_norms_perm' must be an array of np.int32: {_error}", argument="sorted_paralog_norms_perm") from None
    if sorted_paralog_norms_perm.ndim != 1:
        raise _ToxArgumentError(f"'sorted_paralog_norms_perm' must have 1 dimension, but has {sorted_paralog_norms_perm.ndim}", argument="sorted_paralog_norms_perm")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_paralog_analysis_impl import calc_work_arr_paralog_subsets_size
    _calc_work_arr_paralog_subsets_size_result = calc_work_arr_paralog_subsets_size(max_subset_size=max_subset_size, n_genes=n_genes, filtered_paralogs_mask=filtered_paralogs_mask)
    max_subset_size = _calc_work_arr_paralog_subsets_size_result["max_subset_size"]
    n_paralog_subsets = _calc_work_arr_paralog_subsets_size_result["work_array_size"]

    # Fortran cannot check that shared extents agree; this can
    if paralog_norms.shape[0] != n_genes:
        raise _ToxArgumentError(f"'paralog_norms' has {paralog_norms.shape[0]} along axis 0, but 'genes' implies n_genes == {n_genes}", argument="paralog_norms")
    if sorted_paralog_norms_perm.shape[0] != n_genes:
        raise _ToxArgumentError(f"'sorted_paralog_norms_perm' has {sorted_paralog_norms_perm.shape[0]} along axis 0, but 'genes' implies n_genes == {n_genes}", argument="sorted_paralog_norms_perm")
    if genes.shape[0] != n_dims:
        raise _ToxArgumentError(f"'genes' has {genes.shape[0]} along axis 0, but 'ancestor' implies n_dims == {n_dims}", argument="genes")

    # outputs and work arrays, which the caller never sees
    n_results = _ctypes.c_int(0)
    work_arr_paralog_subsets = _np.empty((n_mask_chunks, n_paralog_subsets,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.detect_subfunctionalization_c(
        ancestor,
        genes,
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        _ctypes.byref(_ctypes.c_int(n_mask_chunks)),
        _ctypes.byref(n_results),
        _ctypes.byref(_ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        _ctypes.byref(_ctypes.c_int(n_paralog_subsets)),
        _ctypes.byref(_ctypes.c_double(rdi_threshold)),
        paralog_norms,
        sorted_paralog_norms_perm,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETECT_SUBFUNCTIONALIZATION_ARGUMENTS, _DETECT_SUBFUNCTIONALIZATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    work_arr_paralog_subsets.flags.writeable = False

    return {
        "n_results": n_results.value,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
    }

def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions

    This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    threshold : float
        filter threshold
    n_families : int
        number of families
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes
        The minimum valid value is `(n_genes + 31) / 32`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F'), read-only
        bit mask that will have the indices of genes kept by this pattern set to 1, else 0
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
    Generated from the Fortran procedure `tox_paralog_analysis::filter_paralogs_by_pattern_dosage_effect`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_angles = _np.ascontiguousarray(gene_angles, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_angles' must be an array of np.float64: {_error}", argument="gene_angles") from None
    if gene_angles.ndim != 1:
        raise _ToxArgumentError(f"'gene_angles' must have 1 dimension, but has {gene_angles.ndim}", argument="gene_angles")
    try:
        gene_to_fam = _np.ascontiguousarray(gene_to_fam, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_fam' must be an array of np.int32: {_error}", argument="gene_to_fam") from None
    if gene_to_fam.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}", argument="gene_to_fam")

    # what the inputs already say, rather than asking for it again
    n_genes = gene_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'gene_angles' implies n_genes == {n_genes}", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    masks = _np.empty((n_mask_chunks, n_families,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles,
        _ctypes.byref(_ctypes.c_double(threshold)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        _ctypes.byref(_ctypes.c_int(n_mask_chunks)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENTS, _FILTER_PARALOGS_BY_PATTERN_DOSAGE_EFFECT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    masks.flags.writeable = False

    return masks

def filter_paralogs_by_pattern_subfunctionalization(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks,
):
    r"""Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions

    This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.

    Parameters
    ----------
    gene_angles : np.ndarray[np.float64] of shape (n_genes,)
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
    threshold : float
        filter threshold
    n_families : int
        number of families
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
    n_mask_chunks : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes
        The minimum valid value is `(n_genes + 31) / 32`.

    Returns
    -------
    masks : np.ndarray[np.int32] of shape (n_mask_chunks, n_families,), column-major (order='F'), read-only
        bit mask that will have the indices of genes kept by this pattern set to 1, else 0
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
    Generated from the Fortran procedure `tox_paralog_analysis::filter_paralogs_by_pattern_subfunctionalization`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_angles = _np.ascontiguousarray(gene_angles, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_angles' must be an array of np.float64: {_error}", argument="gene_angles") from None
    if gene_angles.ndim != 1:
        raise _ToxArgumentError(f"'gene_angles' must have 1 dimension, but has {gene_angles.ndim}", argument="gene_angles")
    try:
        gene_to_fam = _np.ascontiguousarray(gene_to_fam, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_fam' must be an array of np.int32: {_error}", argument="gene_to_fam") from None
    if gene_to_fam.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}", argument="gene_to_fam")

    # what the inputs already say, rather than asking for it again
    n_genes = gene_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'gene_angles' implies n_genes == {n_genes}", argument="gene_to_fam")

    # outputs and work arrays, which the caller never sees
    masks = _np.empty((n_mask_chunks, n_families,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.filter_paralogs_by_pattern_subfunctionalization_c(
        gene_angles,
        _ctypes.byref(_ctypes.c_double(threshold)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        _ctypes.byref(_ctypes.c_int(n_mask_chunks)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENTS, _FILTER_PARALOGS_BY_PATTERN_SUBFUNCTIONALIZATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    masks.flags.writeable = False

    return masks
