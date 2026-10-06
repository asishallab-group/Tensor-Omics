r"""tox_get_outliers_by_angle

Gene outliers, from how far each gene's direction deviates from its family's.

Two variants answer the same question on different data. The spherical variant takes
expression vectors: every gene points somewhere in expression space, its family has a mean
direction, and a gene's angular deviation is the angle between the two. The RAP variant takes
one signed angle per gene within a relative axis plane, as
:func:`tensor_omics.clock_hand_angles_for_shift_vectors`
measures them: the family has a circular mean angle, and a gene's angular deviation is its
absolute wrapped distance to that mean.

Both then share the rest. A family's spread is its angular dispersion
\(\sigma = \sqrt{-2 \ln R}\), where \(R\) is the mean resultant length of its members' unit
vectors. A gene's relative angular deviation is its angular deviation divided by its family's
angular dispersion, so a gene in a tight family counts as far off sooner than one in a loose
family. The genes whose relative angular deviation lies strictly above a quantile of all of
them are flagged: a ranking cut, not a statistical test, and values equal to the quantile are
not flagged.

`detect_angle_outliers` and `detect_angle_outliers_rap` run the whole method in one call and
are the entry points to reach for first. Each step is callable on its own too:
`compute_family_direction` (or `compute_family_direction_rap`), `compute_angular_deviations`
(or `compute_angular_deviations_rap`), then the shared `compute_relative_angular_deviations`,
`compute_angle_outlier_threshold` and `flag_angle_outliers`.

A family's statistics exist only when at least three of its genes have a direction and their
unit vectors do not cancel out. A family for which that fails, or whose dispersion lies
outside the accepted bounds, is reported in `status` (without a direction unless its only fault
is too little angular variation) with a `STAT_*` code from
``tox_errors``, and its genes get no relative angular deviation and are never flagged.
A value that does not exist is marked with a sentinel outside the range of real ones: angular
deviations, dispersions and relative angular deviations are never negative and use `-1.0`;
signed angles and family mean angles lie in \((-\pi, \pi]\), and a missing mean angle is
`-4.0`.

A dispersion too small to tell from rounding counts as no angular variation, even with a
minimum of 0: identical or proportional members compute to deviations of rounding size, not
to exactly zero. The limit is about \((m + 1.5\,n + 16) \cdot 2.2 \times 10^{-16}\) rad for
a family of \(m\) members in \(n\) axes (no axis term for signed angles) -- about
\(10^{-14}\) rad for small families, still only about \(2 \times 10^{-11}\) rad for \(10^5\) members.

Python binding, generated from tox_get_outliers_by_angle. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_family_direction_c.restype = None
_lib.compute_family_direction_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_FAMILY_DIRECTION_ARGUMENTS = ("n_axes", "n_genes", "n_families", "expression_vectors", "gene_to_fam", "family_directions", "angular_dispersions", "member_counts", "status", "min_angular_dispersion", "max_angular_dispersion", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_FAMILY_DIRECTION_ARGUMENT_SOURCES = ("expression_vectors", "expression_vectors", "family_directions", None, None, None, None, None, None, None, None, None,)

_lib.compute_angular_deviations_c.restype = None
_lib.compute_angular_deviations_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ANGULAR_DEVIATIONS_ARGUMENTS = ("n_axes", "n_genes", "n_families", "expression_vectors", "family_directions", "gene_to_fam", "angular_deviations", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ANGULAR_DEVIATIONS_ARGUMENT_SOURCES = ("expression_vectors", "expression_vectors", "family_directions", None, None, None, None, None,)

_lib.compute_family_direction_rap_c.restype = None
_lib.compute_family_direction_rap_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_FAMILY_DIRECTION_RAP_ARGUMENTS = ("n_genes", "n_families", "signed_angles", "gene_to_fam", "family_mean_angles", "angular_dispersions", "member_counts", "status", "min_angular_dispersion", "max_angular_dispersion", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_FAMILY_DIRECTION_RAP_ARGUMENT_SOURCES = ("signed_angles", "family_mean_angles", None, None, None, None, None, None, None, None, None,)

_lib.compute_angular_deviations_rap_c.restype = None
_lib.compute_angular_deviations_rap_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ANGULAR_DEVIATIONS_RAP_ARGUMENTS = ("n_genes", "n_families", "signed_angles", "family_mean_angles", "gene_to_fam", "angular_deviations", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ANGULAR_DEVIATIONS_RAP_ARGUMENT_SOURCES = ("signed_angles", "family_mean_angles", None, None, None, None, None,)

_lib.compute_relative_angular_deviations_c.restype = None
_lib.compute_relative_angular_deviations_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RELATIVE_ANGULAR_DEVIATIONS_ARGUMENTS = ("n_genes", "n_families", "angular_deviations", "angular_dispersions", "gene_to_fam", "relative_angular_deviations", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_RELATIVE_ANGULAR_DEVIATIONS_ARGUMENT_SOURCES = ("angular_deviations", "angular_dispersions", None, None, None, None, None,)

_lib.compute_angle_outlier_threshold_c.restype = None
_lib.compute_angle_outlier_threshold_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ANGLE_OUTLIER_THRESHOLD_ARGUMENTS = ("n_genes", "relative_angular_deviations", "threshold", "quantile_level", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ANGLE_OUTLIER_THRESHOLD_ARGUMENT_SOURCES = ("relative_angular_deviations", None, None, None, None,)

_lib.flag_angle_outliers_c.restype = None
_lib.flag_angle_outliers_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FLAG_ANGLE_OUTLIERS_ARGUMENTS = ("n_genes", "relative_angular_deviations", "threshold", "is_outlier", "ierr",)
#: For a derived argument, the one the caller passed it in
_FLAG_ANGLE_OUTLIERS_ARGUMENT_SOURCES = ("relative_angular_deviations", None, None, None, None,)

_lib.detect_angle_outliers_c.restype = None
_lib.detect_angle_outliers_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_ANGLE_OUTLIERS_ARGUMENTS = ("n_axes", "n_genes", "n_families", "expression_vectors", "gene_to_fam", "family_directions", "angular_dispersions", "member_counts", "status", "relative_angular_deviations", "threshold", "is_outlier", "gene_status", "quantile_level", "min_angular_dispersion", "max_angular_dispersion", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETECT_ANGLE_OUTLIERS_ARGUMENT_SOURCES = ("expression_vectors", "expression_vectors", "family_directions", None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.detect_angle_outliers_rap_c.restype = None
_lib.detect_angle_outliers_rap_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETECT_ANGLE_OUTLIERS_RAP_ARGUMENTS = ("n_genes", "n_families", "signed_angles", "gene_to_fam", "family_mean_angles", "angular_dispersions", "member_counts", "status", "relative_angular_deviations", "threshold", "is_outlier", "gene_status", "quantile_level", "min_angular_dispersion", "max_angular_dispersion", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETECT_ANGLE_OUTLIERS_RAP_ARGUMENT_SOURCES = ("signed_angles", "family_mean_angles", None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def compute_family_direction(
        n_families,
        expression_vectors,
        gene_to_fam,
        min_angular_dispersion=0.0,
        max_angular_dispersion=1.1774100225154747,
):
    r"""Compute each family's mean direction and angular dispersion from expression vectors

    Every gene's expression vector is scaled to unit length, and a family's mean direction is
    the normalised sum of its members' unit vectors. With \(R\) the length of that sum divided
    by the number of members, the family's angular dispersion is \(\sigma = \sqrt{-2 \ln R}\):
    zero when all members point the same way, growing as they spread.

    The dispersion is computed from the members' own angular deviations -- the angles
    `compute_angular_deviations` reports -- as \(1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)\),
    which equals the definition because the direction is the members' mean direction, and which
    loses nothing to cancellation when the members are close.

    A gene whose expression vector is zero has no direction. It is left out of its family's
    statistics entirely -- neither counted nor summed. A vector whose components are all
    subnormal (below `tiny(1.0)`, about `2.2e-308`, in magnitude) counts as zero too, in every
    build: whether such numbers survive at all depends on the build's floating-point mode. A
    family whose statistics do not exist gets the zero vector as its direction; see `status`.

    Parameters
    ----------
    n_families : int
        Number of gene families
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        Expression vector of each gene, one per column
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    min_angular_dispersion : float, optional, default 0.0
        Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
        no dispersion is accepted and every family with statistics is reported.
        The default value is `0.0`.
        The minimum valid value is `0.0`.
    max_angular_dispersion : float, optional, default 1.1774100225154747
        Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
        The default value is `sqrt(-2.0*log(0.5))`.
        The minimum valid value is `0.0`.
        The maximum valid value is `5.0`.
        At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

    Returns
    -------
    dict
        with keys:

        family_directions : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
            Mean direction of each family as a unit vector, one per column; the zero vector
            where the family has fewer than three members with a direction or no stable
            direction (see `status`)
            A result is a value; call `.copy()` to obtain a modifiable array.
        angular_dispersions : np.ndarray[np.float64] of shape (n_families,), read-only
            Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            `-1.0` where `status` reports why there is none
            A result is a value; call `.copy()` to obtain a modifiable array.
        member_counts : np.ndarray[np.int32] of shape (n_families,), read-only
            Number of genes of each family that have a direction, i.e. whose expression vector
            is not zero -- the members the statistics are taken over
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_families,), read-only
            Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. ``STAT_TOO_FEW_MEMBERS`` where fewer than three genes of the family have a direction: no direction, no dispersion. ``STAT_NO_STABLE_DIRECTION`` where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. ``STAT_NO_ANGULAR_VARIATION`` where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_family_direction`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    family_directions = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    angular_dispersions = np.empty((n_families,), dtype=np.float64, order='C')
    member_counts = np.empty((n_families,), dtype=np.int32, order='C')
    status = np.empty((n_families,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_family_direction_c(
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        expression_vectors,
        gene_to_fam,
        family_directions,
        angular_dispersions,
        member_counts,
        status,
        ctypes.byref(ctypes.c_double(min_angular_dispersion)),
        ctypes.byref(ctypes.c_double(max_angular_dispersion)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_FAMILY_DIRECTION_ARGUMENTS, _COMPUTE_FAMILY_DIRECTION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    family_directions.flags.writeable = False
    angular_dispersions.flags.writeable = False
    member_counts.flags.writeable = False
    status.flags.writeable = False

    return {
        "family_directions": family_directions,
        "angular_dispersions": angular_dispersions,
        "member_counts": member_counts,
        "status": status,
    }

def compute_angular_deviations(
        expression_vectors,
        family_directions,
        gene_to_fam,
):
    r"""Compute the angle between each gene's expression vector and its family's mean direction

    The angle lies in \([0, \pi]\) and does not depend on the length of either vector.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        Expression vector of each gene, one per column
    family_directions : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        Mean direction of each family, one per column, as
        :func:`tensor_omics.compute_family_direction`
        gives it; the zero vector marks a family without one
        The minimum valid value is `-1.0`.
        The maximum valid value is `1.0`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.

    Returns
    -------
    angular_deviations : np.ndarray[np.float64] of shape (n_genes,), read-only
        Angle in radians between each gene and its family's direction, in \([0, \pi]\), or
        `-1.0` where the gene has no family, its expression vector
        is zero or all subnormal, or its family has no direction
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_angular_deviations`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        family_directions = np.asfortranarray(family_directions, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'family_directions' must be an array of np.float64: {error}") from None
    if family_directions.ndim != 2:
        raise ValueError(f"'family_directions' must have 2 dimensions, but has {family_directions.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_families = family_directions.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if family_directions.shape[0] != n_axes:
        raise ValueError(f"'family_directions' has {family_directions.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_axes == {n_axes}"
        )
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    angular_deviations = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_angular_deviations_c(
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        expression_vectors,
        family_directions,
        gene_to_fam,
        angular_deviations,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ANGULAR_DEVIATIONS_ARGUMENTS, _COMPUTE_ANGULAR_DEVIATIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    angular_deviations.flags.writeable = False

    return angular_deviations

def compute_family_direction_rap(
        n_families,
        signed_angles,
        gene_to_fam,
        min_angular_dispersion=0.0,
        max_angular_dispersion=1.1774100225154747,
):
    r"""Compute each family's circular mean angle and angular dispersion from signed angles

    Every angle is taken as a unit vector \((\cos\theta, \sin\theta)\). A family's mean angle is
    the angle of their sum, and with \(R\) the length of that sum divided by the number of
    members, its angular dispersion is \(\sigma = \sqrt{-2 \ln R}\): zero when all members have
    the same angle, growing as they spread. It is computed from the members' own angular
    deviations -- the distances `compute_angular_deviations_rap` reports -- as
    \(1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)\), which loses nothing to cancellation when
    the members are close.

    Parameters
    ----------
    n_families : int
        Number of gene families
    signed_angles : np.ndarray[np.float64] of shape (n_genes,)
        Signed angle of each gene in radians, in \((-\pi, \pi]\), as
        :func:`tensor_omics.clock_hand_angles_for_shift_vectors`
        measures it. A gene without an angle is one with no family in `gene_to_fam`; its
        entry here is ignored, but must still lie in the range.
        The minimum valid value is `above(-PI)`.
        The maximum valid value is `PI`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    min_angular_dispersion : float, optional, default 0.0
        Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
        no dispersion is accepted and every family with statistics is reported.
        The default value is `0.0`.
        The minimum valid value is `0.0`.
    max_angular_dispersion : float, optional, default 1.1774100225154747
        Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
        The default value is `sqrt(-2.0*log(0.5))`.
        The minimum valid value is `0.0`.
        The maximum valid value is `5.0`.
        At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

    Returns
    -------
    dict
        with keys:

        family_mean_angles : np.ndarray[np.float64] of shape (n_families,), read-only
            Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            `-4.0` where the family has fewer than three members or no
            stable direction (see `status`)
            A result is a value; call `.copy()` to obtain a modifiable array.
        angular_dispersions : np.ndarray[np.float64] of shape (n_families,), read-only
            Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            `-1.0` where `status` reports why there is none
            A result is a value; call `.copy()` to obtain a modifiable array.
        member_counts : np.ndarray[np.int32] of shape (n_families,), read-only
            Number of genes of each family -- the members the statistics are taken over
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_families,), read-only
            Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. ``STAT_TOO_FEW_MEMBERS`` where fewer than three genes belong to the family: no mean angle, no dispersion. ``STAT_NO_STABLE_DIRECTION`` where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. ``STAT_NO_ANGULAR_VARIATION`` where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_family_direction_rap`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        signed_angles = np.ascontiguousarray(signed_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'signed_angles' must be an array of np.float64: {error}") from None
    if signed_angles.ndim != 1:
        raise ValueError(f"'signed_angles' must have 1 dimension, but has {signed_angles.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = signed_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'signed_angles' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    family_mean_angles = np.empty((n_families,), dtype=np.float64, order='C')
    angular_dispersions = np.empty((n_families,), dtype=np.float64, order='C')
    member_counts = np.empty((n_families,), dtype=np.int32, order='C')
    status = np.empty((n_families,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_family_direction_rap_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        signed_angles,
        gene_to_fam,
        family_mean_angles,
        angular_dispersions,
        member_counts,
        status,
        ctypes.byref(ctypes.c_double(min_angular_dispersion)),
        ctypes.byref(ctypes.c_double(max_angular_dispersion)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_FAMILY_DIRECTION_RAP_ARGUMENTS, _COMPUTE_FAMILY_DIRECTION_RAP_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    family_mean_angles.flags.writeable = False
    angular_dispersions.flags.writeable = False
    member_counts.flags.writeable = False
    status.flags.writeable = False

    return {
        "family_mean_angles": family_mean_angles,
        "angular_dispersions": angular_dispersions,
        "member_counts": member_counts,
        "status": status,
    }

def compute_angular_deviations_rap(
        signed_angles,
        family_mean_angles,
        gene_to_fam,
):
    r"""Compute the absolute angular distance between each gene's signed angle and its family's mean angle

    The distance is taken around the circle, so it lies in \([0, \pi]\): an angle just below
    \(\pi\) and one just above \(-\pi\) are close. See ``wrap_angle``.

    Parameters
    ----------
    signed_angles : np.ndarray[np.float64] of shape (n_genes,)
        Signed angle of each gene in radians, in \((-\pi, \pi]\), as
        :func:`tensor_omics.clock_hand_angles_for_shift_vectors`
        measures it. A gene without an angle is one with no family in `gene_to_fam`; its
        entry here is ignored, but must still lie in the range.
        The minimum valid value is `above(-PI)`.
        The maximum valid value is `PI`.
    family_mean_angles : np.ndarray[np.float64] of shape (n_families,)
        Circular mean angle of each family in radians, as
        :func:`tensor_omics.compute_family_direction_rap`
        gives it, or `-4.0` for a family without one
        The minimum valid value is `above(-PI)`.
        The maximum valid value is `PI`.
        The value `-4.0` is additionally accepted.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.

    Returns
    -------
    angular_deviations : np.ndarray[np.float64] of shape (n_genes,), read-only
        Absolute angular distance in radians between each gene's angle and its family's mean
        angle, in \([0, \pi]\), or `-1.0` where the gene has no
        family or its family has no mean angle
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_angular_deviations_rap`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        signed_angles = np.ascontiguousarray(signed_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'signed_angles' must be an array of np.float64: {error}") from None
    if signed_angles.ndim != 1:
        raise ValueError(f"'signed_angles' must have 1 dimension, but has {signed_angles.ndim}")
    try:
        family_mean_angles = np.ascontiguousarray(family_mean_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'family_mean_angles' must be an array of np.float64: {error}") from None
    if family_mean_angles.ndim != 1:
        raise ValueError(f"'family_mean_angles' must have 1 dimension, but has {family_mean_angles.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = signed_angles.shape[0]
    n_families = family_mean_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'signed_angles' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    angular_deviations = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_angular_deviations_rap_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        signed_angles,
        family_mean_angles,
        gene_to_fam,
        angular_deviations,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ANGULAR_DEVIATIONS_RAP_ARGUMENTS, _COMPUTE_ANGULAR_DEVIATIONS_RAP_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    angular_deviations.flags.writeable = False

    return angular_deviations

def compute_relative_angular_deviations(
        angular_deviations,
        angular_dispersions,
        gene_to_fam,
):
    r"""Divide each gene's angular deviation by its family's angular dispersion

    The relative angular deviation says how far off a gene is in units of its own family's
    spread, so genes of tight and loose families can be ranked together. It is a plain quotient,
    with nothing subtracted. It serves both variants: the angular deviations may come from
    :func:`tensor_omics.compute_angular_deviations` or
    from :func:`tensor_omics.compute_angular_deviations_rap`.

    Parameters
    ----------
    angular_deviations : np.ndarray[np.float64] of shape (n_genes,)
        Angular deviation of each gene from its family in radians, or
        `-1.0` where there is none
        The minimum valid value is `0.0`.
        The maximum valid value is `PI`.
        The value `-1.0` is additionally accepted.
    angular_dispersions : np.ndarray[np.float64] of shape (n_families,)
        Angular dispersion of each family, or `-1.0` where there
        is none
        The minimum valid value is `0.0`.
        The value `-1.0` is additionally accepted.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `angular_deviations`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.

    Returns
    -------
    relative_angular_deviations : np.ndarray[np.float64] of shape (n_genes,), read-only
        Angular deviation divided by the family's angular dispersion, at least zero, or
        `-1.0` where the gene has no family or no
        angular deviation, or its family's dispersion is the sentinel, zero, or below the
        smallest normal number `tiny(1.0)` (about `2.2e-308`) -- a dispersion so small that
        the quotient could exceed the largest representable number
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_relative_angular_deviations`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        angular_deviations = np.ascontiguousarray(angular_deviations, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'angular_deviations' must be an array of np.float64: {error}") from None
    if angular_deviations.ndim != 1:
        raise ValueError(f"'angular_deviations' must have 1 dimension, but has {angular_deviations.ndim}")
    try:
        angular_dispersions = np.ascontiguousarray(angular_dispersions, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'angular_dispersions' must be an array of np.float64: {error}") from None
    if angular_dispersions.ndim != 1:
        raise ValueError(f"'angular_dispersions' must have 1 dimension, but has {angular_dispersions.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = angular_deviations.shape[0]
    n_families = angular_dispersions.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'angular_deviations' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    relative_angular_deviations = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_relative_angular_deviations_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        angular_deviations,
        angular_dispersions,
        gene_to_fam,
        relative_angular_deviations,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_RELATIVE_ANGULAR_DEVIATIONS_ARGUMENTS, _COMPUTE_RELATIVE_ANGULAR_DEVIATIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    relative_angular_deviations.flags.writeable = False

    return relative_angular_deviations

def compute_angle_outlier_threshold(
        relative_angular_deviations,
        quantile_level=0.95,
):
    r"""Compute the threshold at which a relative angular deviation counts as an outlier

    The threshold is the `quantile_level` quantile of the relative angular deviations that
    exist, interpolated linearly between neighbouring values (the "type 7" quantile). The
    sentinels do not take part.

    This is a ranking cut, not a statistical test: whatever the data, it picks out the largest
    values. :func:`tensor_omics.flag_angle_outliers`
    flags only values strictly above the threshold, so with the default level 0.95 at most about
    the top 5% are flagged. Where many values tie at the quantile position -- a family of
    identical genes plus one outlier -- the threshold equals the tied value and only the values
    above it are flagged; where all values are equal, nothing is. With fewer than two values
    there is nothing to rank, and the threshold is the largest representable number, so nothing
    is flagged.

    Parameters
    ----------
    relative_angular_deviations : np.ndarray[np.float64] of shape (n_genes,)
        Relative angular deviation of each gene, or `-1.0`
        where there is none
        The minimum valid value is `0.0`.
        The value `-1.0` is additionally accepted.
    quantile_level : float, optional, default 0.95
        Quantile level in \([0, 1]\) the threshold is taken at
        The default value is `0.95`.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    threshold : float
        The relative angular deviation a gene must exceed to be an outlier

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::compute_angle_outlier_threshold`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        relative_angular_deviations = np.ascontiguousarray(relative_angular_deviations, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'relative_angular_deviations' must be an array of np.float64: {error}") from None
    if relative_angular_deviations.ndim != 1:
        raise ValueError(f"'relative_angular_deviations' must have 1 dimension, but has {relative_angular_deviations.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = relative_angular_deviations.shape[0]

    # outputs and work arrays, which the caller never sees
    threshold = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.compute_angle_outlier_threshold_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        relative_angular_deviations,
        ctypes.byref(threshold),
        ctypes.byref(ctypes.c_double(quantile_level)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ANGLE_OUTLIER_THRESHOLD_ARGUMENTS, _COMPUTE_ANGLE_OUTLIER_THRESHOLD_ARGUMENT_SOURCES)

    return threshold.value

def flag_angle_outliers(
        relative_angular_deviations,
        threshold,
):
    r"""Flag the genes whose relative angular deviation lies above the threshold

    A gene is an outlier when its relative angular deviation is positive and strictly greater
    than `threshold`; a value equal to the threshold is not flagged. So a gene without a value,
    or with a deviation of zero, is never flagged, whatever the threshold.

    Parameters
    ----------
    relative_angular_deviations : np.ndarray[np.float64] of shape (n_genes,)
        Relative angular deviation of each gene, or `-1.0`
        where there is none
        The minimum valid value is `0.0`.
        The value `-1.0` is additionally accepted.
    threshold : float
        The relative angular deviation a gene must exceed to be an outlier, as
        :func:`tensor_omics.compute_angle_outlier_threshold`
        computes it

    Returns
    -------
    is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
        `True` for each gene that is an outlier
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::flag_angle_outliers`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        relative_angular_deviations = np.ascontiguousarray(relative_angular_deviations, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'relative_angular_deviations' must be an array of np.float64: {error}") from None
    if relative_angular_deviations.ndim != 1:
        raise ValueError(f"'relative_angular_deviations' must have 1 dimension, but has {relative_angular_deviations.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = relative_angular_deviations.shape[0]

    # outputs and work arrays, which the caller never sees
    is_outlier = np.empty((n_genes,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.flag_angle_outliers_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        relative_angular_deviations,
        ctypes.byref(ctypes.c_double(threshold)),
        is_outlier,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FLAG_ANGLE_OUTLIERS_ARGUMENTS, _FLAG_ANGLE_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    is_outlier.flags.writeable = False

    return is_outlier

def detect_angle_outliers(
        n_families,
        expression_vectors,
        gene_to_fam,
        quantile_level=0.95,
        min_angular_dispersion=0.0,
        max_angular_dispersion=1.1774100225154747,
):
    r"""Detect angle outliers among genes from their expression vectors

    Runs the spherical variant end to end:
    :func:`tensor_omics.compute_family_direction`,
    :func:`tensor_omics.compute_angular_deviations`,
    :func:`tensor_omics.compute_relative_angular_deviations`,
    :func:`tensor_omics.compute_angle_outlier_threshold` and
    :func:`tensor_omics.flag_angle_outliers`.

    The threshold is a ranking cut, not a statistical test, and only values strictly above it are
    flagged: where genes tie at the quantile position (identical genes plus one outlier), only
    the values above the tie are flagged, and where all are equal, none is.

    Parameters
    ----------
    n_families : int
        Number of gene families
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        Expression vector of each gene, one per column
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    quantile_level : float, optional, default 0.95
        Quantile level in \([0, 1]\) the threshold is taken at
        The default value is `0.95`.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    min_angular_dispersion : float, optional, default 0.0
        Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
        no dispersion is accepted and every family with statistics is reported.
        The default value is `0.0`.
        The minimum valid value is `0.0`.
    max_angular_dispersion : float, optional, default 1.1774100225154747
        Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
        The default value is `sqrt(-2.0*log(0.5))`.
        The minimum valid value is `0.0`.
        The maximum valid value is `5.0`.
        At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

    Returns
    -------
    dict
        with keys:

        family_directions : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
            Mean direction of each family as a unit vector, one per column; the zero vector
            where the family has fewer than three members with a direction or no stable
            direction (see `status`)
            A result is a value; call `.copy()` to obtain a modifiable array.
        angular_dispersions : np.ndarray[np.float64] of shape (n_families,), read-only
            Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            `-1.0` where `status` reports why there is none
            A result is a value; call `.copy()` to obtain a modifiable array.
        member_counts : np.ndarray[np.int32] of shape (n_families,), read-only
            Number of genes of each family that have a direction, i.e. whose expression vector
            is not zero -- the members the statistics are taken over
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_families,), read-only
            Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. ``STAT_TOO_FEW_MEMBERS`` where fewer than three genes of the family have a direction: no direction, no dispersion. ``STAT_NO_STABLE_DIRECTION`` where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. ``STAT_NO_ANGULAR_VARIATION`` where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.
            A result is a value; call `.copy()` to obtain a modifiable array.
        relative_angular_deviations : np.ndarray[np.float64] of shape (n_genes,), read-only
            Angle between each gene and its family's direction, divided by the family's angular
            dispersion; `-1.0` where `gene_status` is not zero
            A result is a value; call `.copy()` to obtain a modifiable array.
        threshold : float
            The relative angular deviation a gene must exceed to be an outlier; the largest
            representable number when fewer than two genes have a relative angular deviation
        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            `True` for each gene whose relative angular deviation is positive and strictly above
            `threshold`
            A result is a value; call `.copy()` to obtain a modifiable array.
        gene_status : np.ndarray[np.int32] of shape (n_genes,), read-only
            Zero where the gene has a relative angular deviation, otherwise the reason it has
            not: ``STAT_NO_FAMILY`` for a gene without a family,
            ``STAT_ZERO_VECTOR`` for one whose expression vector is
            zero or all subnormal, and otherwise its family's `status`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::detect_angle_outliers`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    family_directions = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    angular_dispersions = np.empty((n_families,), dtype=np.float64, order='C')
    member_counts = np.empty((n_families,), dtype=np.int32, order='C')
    status = np.empty((n_families,), dtype=np.int32, order='C')
    relative_angular_deviations = np.empty((n_genes,), dtype=np.float64, order='C')
    threshold = ctypes.c_double(0)
    is_outlier = np.empty((n_genes,), dtype=np.bool_, order='C')
    gene_status = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.detect_angle_outliers_c(
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        expression_vectors,
        gene_to_fam,
        family_directions,
        angular_dispersions,
        member_counts,
        status,
        relative_angular_deviations,
        ctypes.byref(threshold),
        is_outlier,
        gene_status,
        ctypes.byref(ctypes.c_double(quantile_level)),
        ctypes.byref(ctypes.c_double(min_angular_dispersion)),
        ctypes.byref(ctypes.c_double(max_angular_dispersion)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_ANGLE_OUTLIERS_ARGUMENTS, _DETECT_ANGLE_OUTLIERS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    family_directions.flags.writeable = False
    angular_dispersions.flags.writeable = False
    member_counts.flags.writeable = False
    status.flags.writeable = False
    relative_angular_deviations.flags.writeable = False
    is_outlier.flags.writeable = False
    gene_status.flags.writeable = False

    return {
        "family_directions": family_directions,
        "angular_dispersions": angular_dispersions,
        "member_counts": member_counts,
        "status": status,
        "relative_angular_deviations": relative_angular_deviations,
        "threshold": threshold.value,
        "is_outlier": is_outlier,
        "gene_status": gene_status,
    }

def detect_angle_outliers_rap(
        n_families,
        signed_angles,
        gene_to_fam,
        quantile_level=0.95,
        min_angular_dispersion=0.0,
        max_angular_dispersion=1.1774100225154747,
):
    r"""Detect angle outliers among genes from their signed angles in a relative axis plane

    Runs the RAP variant end to end:
    :func:`tensor_omics.compute_family_direction_rap`,
    :func:`tensor_omics.compute_angular_deviations_rap`,
    :func:`tensor_omics.compute_relative_angular_deviations`,
    :func:`tensor_omics.compute_angle_outlier_threshold` and
    :func:`tensor_omics.flag_angle_outliers`.

    The threshold is a ranking cut, not a statistical test, and only values strictly above it are
    flagged: where genes tie at the quantile position (identical genes plus one outlier), only
    the values above the tie are flagged, and where all are equal, none is.

    Parameters
    ----------
    n_families : int
        Number of gene families
    signed_angles : np.ndarray[np.float64] of shape (n_genes,)
        Signed angle of each gene in radians, in \((-\pi, \pi]\), as
        :func:`tensor_omics.clock_hand_angles_for_shift_vectors`
        measures it. A gene without an angle is one with no family in `gene_to_fam`; its
        entry here is ignored, but must still lie in the range.
        The minimum valid value is `above(-PI)`.
        The maximum valid value is `PI`.
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    quantile_level : float, optional, default 0.95
        Quantile level in \([0, 1]\) the threshold is taken at
        The default value is `0.95`.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    min_angular_dispersion : float, optional, default 0.0
        Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
        no dispersion is accepted and every family with statistics is reported.
        The default value is `0.0`.
        The minimum valid value is `0.0`.
    max_angular_dispersion : float, optional, default 1.1774100225154747
        Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
        The default value is `sqrt(-2.0*log(0.5))`.
        The minimum valid value is `0.0`.
        The maximum valid value is `5.0`.
        At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

    Returns
    -------
    dict
        with keys:

        family_mean_angles : np.ndarray[np.float64] of shape (n_families,), read-only
            Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            `-4.0` where the family has fewer than three members or no
            stable direction (see `status`)
            A result is a value; call `.copy()` to obtain a modifiable array.
        angular_dispersions : np.ndarray[np.float64] of shape (n_families,), read-only
            Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            `-1.0` where `status` reports why there is none
            A result is a value; call `.copy()` to obtain a modifiable array.
        member_counts : np.ndarray[np.int32] of shape (n_families,), read-only
            Number of genes of each family -- the members the statistics are taken over
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_families,), read-only
            Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. ``STAT_TOO_FEW_MEMBERS`` where fewer than three genes belong to the family: no mean angle, no dispersion. ``STAT_NO_STABLE_DIRECTION`` where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. ``STAT_NO_ANGULAR_VARIATION`` where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.
            A result is a value; call `.copy()` to obtain a modifiable array.
        relative_angular_deviations : np.ndarray[np.float64] of shape (n_genes,), read-only
            Absolute angular distance between each gene's angle and its family's mean angle,
            divided by the family's angular dispersion; `-1.0`
            where `gene_status` is not zero
            A result is a value; call `.copy()` to obtain a modifiable array.
        threshold : float
            The relative angular deviation a gene must exceed to be an outlier; the largest
            representable number when fewer than two genes have a relative angular deviation
        is_outlier : np.ndarray[np.bool_] of shape (n_genes,), read-only
            `True` for each gene whose relative angular deviation is positive and strictly above
            `threshold`
            A result is a value; call `.copy()` to obtain a modifiable array.
        gene_status : np.ndarray[np.int32] of shape (n_genes,), read-only
            Zero where the gene has a relative angular deviation, otherwise the reason it has
            not: ``STAT_NO_FAMILY`` for a gene without a family, and
            otherwise its family's `status`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_get_outliers_by_angle::detect_angle_outliers_rap`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        signed_angles = np.ascontiguousarray(signed_angles, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'signed_angles' must be an array of np.float64: {error}") from None
    if signed_angles.ndim != 1:
        raise ValueError(f"'signed_angles' must have 1 dimension, but has {signed_angles.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = signed_angles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'signed_angles' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    family_mean_angles = np.empty((n_families,), dtype=np.float64, order='C')
    angular_dispersions = np.empty((n_families,), dtype=np.float64, order='C')
    member_counts = np.empty((n_families,), dtype=np.int32, order='C')
    status = np.empty((n_families,), dtype=np.int32, order='C')
    relative_angular_deviations = np.empty((n_genes,), dtype=np.float64, order='C')
    threshold = ctypes.c_double(0)
    is_outlier = np.empty((n_genes,), dtype=np.bool_, order='C')
    gene_status = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.detect_angle_outliers_rap_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        signed_angles,
        gene_to_fam,
        family_mean_angles,
        angular_dispersions,
        member_counts,
        status,
        relative_angular_deviations,
        ctypes.byref(threshold),
        is_outlier,
        gene_status,
        ctypes.byref(ctypes.c_double(quantile_level)),
        ctypes.byref(ctypes.c_double(min_angular_dispersion)),
        ctypes.byref(ctypes.c_double(max_angular_dispersion)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETECT_ANGLE_OUTLIERS_RAP_ARGUMENTS, _DETECT_ANGLE_OUTLIERS_RAP_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    family_mean_angles.flags.writeable = False
    angular_dispersions.flags.writeable = False
    member_counts.flags.writeable = False
    status.flags.writeable = False
    relative_angular_deviations.flags.writeable = False
    is_outlier.flags.writeable = False
    gene_status.flags.writeable = False

    return {
        "family_mean_angles": family_mean_angles,
        "angular_dispersions": angular_dispersions,
        "member_counts": member_counts,
        "status": status,
        "relative_angular_deviations": relative_angular_deviations,
        "threshold": threshold.value,
        "is_outlier": is_outlier,
        "gene_status": gene_status,
    }
