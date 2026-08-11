"""tox_flyer_json

Serialization of tox analysis results into the JSON format consumed by the tox_flyer viewer.

Python binding, generated from tox_flyer_json. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.serialize_tox_data_as_flyer_json_c.restype = None
_lib.serialize_tox_data_as_flyer_json_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_TOX_DATA_AS_FLYER_JSON_ARGUMENTS = ("filename", "tissues", "n_tissues", "family_ids", "n_families", "centroids", "gene_ids", "n_genes", "genes", "gene_to_fam", "sorted_gene_to_fam_perm", "gene_outliers", "gene_species", "gene_types", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_TOX_DATA_AS_FLYER_JSON_ARGUMENT_SOURCES = (None, None, "tissues", None, "family_ids", None, None, "gene_ids", None, None, None, None, None, None, None,)

def serialize_tox_data_as_flyer_json(
        filename,
        tissues,
        family_ids,
        centroids,
        gene_ids,
        genes,
        gene_to_fam,
        sorted_gene_to_fam_perm,
        gene_outliers,
        gene_species,
        gene_types,
):
    r"""Serializes tox related data to JSON, compatible with the tox_flyer

    Parameters
    ----------
    filename : str
        Name of the file to write the output to
    tissues : sequence of str, of length n_tissues
        Tissue identifiers
    family_ids : sequence of str, of length n_families
        Family identifiers
    centroids : np.ndarray[np.float64] of shape (n_tissues, n_families,), column-major (order='F')
        Centroid data
    gene_ids : sequence of str, of length n_genes
        Gene identifiers
    genes : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene data
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Gene index to Family index mapping
    sorted_gene_to_fam_perm : np.ndarray[np.int32] of shape (n_genes,)
        Permutation vector that sorts `gene_to_fam`
    gene_outliers : np.ndarray[np.bool_] of shape (n_genes,)
        Specifies if a gene is an outlier
    gene_species : sequence of str, of length n_genes
        Species name per gene
    gene_types : sequence of str, of length n_genes
        Gene type string (ortholog/paralog)

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_flyer_json::serialize_tox_data_as_flyer_json`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        tissues = np.asarray([str(_s).encode() for _s in tissues], dtype="S")
    except TypeError as error:
        raise TypeError(f"'tissues' must be a sequence of strings: {error}") from None
    if tissues.ndim != 1:
        raise ValueError(f"'tissues' must have 1 dimension, but has {tissues.ndim}")
    try:
        family_ids = np.asarray([str(_s).encode() for _s in family_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'family_ids' must be a sequence of strings: {error}") from None
    if family_ids.ndim != 1:
        raise ValueError(f"'family_ids' must have 1 dimension, but has {family_ids.ndim}")
    try:
        centroids = np.asfortranarray(centroids, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'centroids' must be an array of np.float64: {error}") from None
    if centroids.ndim != 2:
        raise ValueError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}")
    try:
        gene_ids = np.asarray([str(_s).encode() for _s in gene_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")
    try:
        genes = np.asfortranarray(genes, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'genes' must be an array of np.float64: {error}") from None
    if genes.ndim != 2:
        raise ValueError(f"'genes' must have 2 dimensions, but has {genes.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
    try:
        sorted_gene_to_fam_perm = np.ascontiguousarray(sorted_gene_to_fam_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_gene_to_fam_perm' must be an array of np.int32: {error}") from None
    if sorted_gene_to_fam_perm.ndim != 1:
        raise ValueError(f"'sorted_gene_to_fam_perm' must have 1 dimension, but has {sorted_gene_to_fam_perm.ndim}")
    try:
        gene_outliers = np.ascontiguousarray(gene_outliers, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_outliers' must be an array of np.bool_: {error}") from None
    if gene_outliers.ndim != 1:
        raise ValueError(f"'gene_outliers' must have 1 dimension, but has {gene_outliers.ndim}")
    try:
        gene_species = np.asarray([str(_s).encode() for _s in gene_species], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_species' must be a sequence of strings: {error}") from None
    if gene_species.ndim != 1:
        raise ValueError(f"'gene_species' must have 1 dimension, but has {gene_species.ndim}")
    try:
        gene_types = np.asarray([str(_s).encode() for _s in gene_types], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_types' must be a sequence of strings: {error}") from None
    if gene_types.ndim != 1:
        raise ValueError(f"'gene_types' must have 1 dimension, but has {gene_types.ndim}")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    tissues_strlen = tissues.itemsize
    n_tissues = tissues.shape[0]
    family_ids_strlen = family_ids.itemsize
    n_families = family_ids.shape[0]
    gene_ids_strlen = gene_ids.itemsize
    n_genes = gene_ids.shape[0]
    gene_species_strlen = gene_species.itemsize
    gene_types_strlen = gene_types.itemsize

    # Fortran cannot check that shared extents agree; this can
    if centroids.shape[0] != n_tissues:
        raise ValueError(f"'centroids' has {centroids.shape[0]} along axis 0, but "
            f"'tissues' implies n_tissues == {n_tissues}"
        )
    if genes.shape[0] != n_tissues:
        raise ValueError(f"'genes' has {genes.shape[0]} along axis 0, but "
            f"'tissues' implies n_tissues == {n_tissues}"
        )
    if centroids.shape[1] != n_families:
        raise ValueError(f"'centroids' has {centroids.shape[1]} along axis 1, but "
            f"'family_ids' implies n_families == {n_families}"
        )
    if genes.shape[1] != n_genes:
        raise ValueError(f"'genes' has {genes.shape[1]} along axis 1, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )
    if sorted_gene_to_fam_perm.shape[0] != n_genes:
        raise ValueError(f"'sorted_gene_to_fam_perm' has {sorted_gene_to_fam_perm.shape[0]} along axis 0, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )
    if gene_outliers.shape[0] != n_genes:
        raise ValueError(f"'gene_outliers' has {gene_outliers.shape[0]} along axis 0, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )
    if gene_species.shape[0] != n_genes:
        raise ValueError(f"'gene_species' has {gene_species.shape[0]} along axis 0, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )
    if gene_types.shape[0] != n_genes:
        raise ValueError(f"'gene_types' has {gene_types.shape[0]} along axis 0, but "
            f"'gene_ids' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_tox_data_as_flyer_json_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        tissues,
        ctypes.byref(ctypes.c_int(tissues_strlen)),
        ctypes.byref(ctypes.c_int(n_tissues)),
        family_ids,
        ctypes.byref(ctypes.c_int(family_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_families)),
        centroids,
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_genes)),
        genes,
        gene_to_fam,
        sorted_gene_to_fam_perm,
        gene_outliers,
        gene_species,
        ctypes.byref(ctypes.c_int(gene_species_strlen)),
        gene_types,
        ctypes.byref(ctypes.c_int(gene_types_strlen)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_TOX_DATA_AS_FLYER_JSON_ARGUMENTS, _SERIALIZE_TOX_DATA_AS_FLYER_JSON_ARGUMENT_SOURCES)

    return None
