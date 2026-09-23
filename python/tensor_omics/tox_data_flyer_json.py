r"""tox_data_flyer_json

Export of a data set to the flyer: a JSON file describing genes, their families and their
expression vectors, for viewing them in three dimensions.

:func:`tensor_omics.save_flyer_json` is the entry point; its
documentation is the specification of the file format.

Python binding, generated from tox_data_flyer_json. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.save_flyer_json_c.restype = None
_lib.save_flyer_json_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SAVE_FLYER_JSON_ARGUMENTS = ("n_axes", "n_genes", "n_families", "expression_vectors", "family_centroids", "gene_to_fam", "is_outlier", "axis_labels", "family_ids", "gene_ids", "gene_species", "gene_types", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_SAVE_FLYER_JSON_ARGUMENT_SOURCES = ("expression_vectors", "expression_vectors", "family_centroids", None, None, None, None, None, None, None, None, None, None, None,)

def save_flyer_json(
        expression_vectors,
        family_centroids,
        gene_to_fam,
        is_outlier,
        axis_labels,
        family_ids,
        gene_ids,
        gene_species,
        gene_types,
        filename,
):
    r"""Write a data set as a flyer JSON file

    Writes the genes, their expression vectors, their families and the family centroids into
    a new JSON file, the input of the flyer. Everything is checked before the file is
    created, and nothing is left behind when a check or a write fails.

    ### The file

    One UTF-8 JSON object on a single line, without whitespace, ended by a line feed. With
    two axes, one family and one gene it reads (broken into lines here):

    ```
    {"kind":"tox-flyer","version":1,"tissues":["liver","brain"],
    "families":[{"family":"F1","gene_indices":[1],"centroid":[5.0000000000000000E-001,-2.0000000000000000E+000]}],
    "genes":[{"coordinates":[1.0000000000000000E+000,0.0],"id":"g1","family":"F1",
    "species":"human","is_outlier":false,"type":"ortholog"}]}
    ```

    The keys of the object:

    - `"kind"`: always `"tox-flyer"`;
    - `"version"`: always `1`;
    - `"tissues"`: `axis_labels` in order, the name of each coordinate axis;
    - `"families"`: one object per family, in the order of `family_ids`;
    - `"genes"`: one object per gene, in input order.

    The keys of each family:

    - `"family"`: its id;
    - `"gene_indices"`: the 1-based positions in `"genes"` of its members, ascending (`[]`
    for a family without genes);
    - `"centroid"`: its column of `family_centroids`, one number per axis.

    The keys of each gene, all six always present:

    - `"coordinates"`: its column of `expression_vectors`, one number per axis;
    - `"id"`: its id;
    - `"family"`: the id of its family, or `null` for a gene in no family (never left out);
    - `"species"`: its species;
    - `"is_outlier"`: `true` or `false`;
    - `"type"`: its type, e.g. ortholog or paralog.

    Readers ignore keys they do not know, so a key can be added without changing the
    version. Ids and labels are unique and non-empty. Numbers carry 17 significant digits,
    which read back as exactly the same double; a zero is written `0.0` or `-0.0`. Strings
    lose their trailing blanks (so a value that genuinely ends in one cannot be written);
    `"`, the backslash and the control characters are escaped.

    ### Empty inputs

    Zero genes give `"genes":[]` and families whose `"gene_indices"` are all `[]`. Zero
    families give `"families":[]`, with every gene's `"family"` null. The Python and R
    bindings cannot pass an empty array of strings yet (issue #221), so there it takes at
    least one family and one gene. At least one axis is always required.

    ### Errors

    Each is reported with the argument it concerns, and the file is not created:

    - `ERR_EMPTY_INPUT`: no axes;
    - `ERR_NAN_INF`: a NaN or an infinity in `expression_vectors` or `family_centroids`;
    - `ERR_INVALID_INPUT`: a `gene_to_fam` entry that is negative or above the number of
    families; an empty or a repeated entry in `axis_labels`, `family_ids` or `gene_ids`;
    - `ERR_INVALID_UTF8`: a string that is not valid UTF-8;
    - `ERR_FILE_OPEN`: `filename` exists already (it is never overwritten) or cannot be
    created.

    `ERR_WRITE_DATA` means writing failed, for example on a full disk; the partial file is
    deleted.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        Coordinates of each gene, one column per gene
    family_centroids : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        Centroid of each family, one column per family
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    is_outlier : np.ndarray[np.bool_] of shape (n_genes,)
        Whether each gene is an outlier
    axis_labels : sequence of str, of length n_axes
        Name of each coordinate axis (e.g. tissues); unique and non-empty
    family_ids : sequence of str, of length n_families
        Id of each family; unique and non-empty
    gene_ids : sequence of str, of length n_genes
        Id of each gene; unique and non-empty
    gene_species : sequence of str, of length n_genes
        Species of each gene
    gene_types : sequence of str, of length n_genes
        Type of each gene, e.g. ortholog or paralog
    filename : str
        Path of the JSON file to create; it must not exist yet

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_flyer_json::save_flyer_json`, whose argument names are
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
        family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'family_centroids' must be an array of np.float64: {error}") from None
    if family_centroids.ndim != 2:
        raise ValueError(f"'family_centroids' must have 2 dimensions, but has {family_centroids.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
    try:
        is_outlier = np.ascontiguousarray(is_outlier, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'is_outlier' must be an array of np.bool_: {error}") from None
    if is_outlier.ndim != 1:
        raise ValueError(f"'is_outlier' must have 1 dimension, but has {is_outlier.ndim}")
    try:
        _axis_labels_bytes = [str(_s).encode() for _s in axis_labels]
        _axis_labels_width = max(map(len, _axis_labels_bytes), default=0) or 1
        axis_labels = np.asarray([_b.ljust(_axis_labels_width) for _b in _axis_labels_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'axis_labels' must be a sequence of strings: {error}") from None
    if axis_labels.ndim != 1:
        raise ValueError(f"'axis_labels' must have 1 dimension, but has {axis_labels.ndim}")
    try:
        _family_ids_bytes = [str(_s).encode() for _s in family_ids]
        _family_ids_width = max(map(len, _family_ids_bytes), default=0) or 1
        family_ids = np.asarray([_b.ljust(_family_ids_width) for _b in _family_ids_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'family_ids' must be a sequence of strings: {error}") from None
    if family_ids.ndim != 1:
        raise ValueError(f"'family_ids' must have 1 dimension, but has {family_ids.ndim}")
    try:
        _gene_ids_bytes = [str(_s).encode() for _s in gene_ids]
        _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
        gene_ids = np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")
    try:
        _gene_species_bytes = [str(_s).encode() for _s in gene_species]
        _gene_species_width = max(map(len, _gene_species_bytes), default=0) or 1
        gene_species = np.asarray([_b.ljust(_gene_species_width) for _b in _gene_species_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_species' must be a sequence of strings: {error}") from None
    if gene_species.ndim != 1:
        raise ValueError(f"'gene_species' must have 1 dimension, but has {gene_species.ndim}")
    try:
        _gene_types_bytes = [str(_s).encode() for _s in gene_types]
        _gene_types_width = max(map(len, _gene_types_bytes), default=0) or 1
        gene_types = np.asarray([_b.ljust(_gene_types_width) for _b in _gene_types_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_types' must be a sequence of strings: {error}") from None
    if gene_types.ndim != 1:
        raise ValueError(f"'gene_types' must have 1 dimension, but has {gene_types.ndim}")
    filename = np.array([str(filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_families = family_centroids.shape[1]
    axis_labels_strlen = axis_labels.itemsize
    family_ids_strlen = family_ids.itemsize
    gene_ids_strlen = gene_ids.itemsize
    gene_species_strlen = gene_species.itemsize
    gene_types_strlen = gene_types.itemsize
    filename_strlen = filename.itemsize

    # Fortran cannot check that shared extents agree; this can
    if family_centroids.shape[0] != n_axes:
        raise ValueError(f"'family_centroids' has {family_centroids.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_axes == {n_axes}"
        )
    if axis_labels.shape[0] != n_axes:
        raise ValueError(f"'axis_labels' has {axis_labels.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_axes == {n_axes}"
        )
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if is_outlier.shape[0] != n_genes:
        raise ValueError(f"'is_outlier' has {is_outlier.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if gene_ids.shape[0] != n_genes:
        raise ValueError(f"'gene_ids' has {gene_ids.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if gene_species.shape[0] != n_genes:
        raise ValueError(f"'gene_species' has {gene_species.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if gene_types.shape[0] != n_genes:
        raise ValueError(f"'gene_types' has {gene_types.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if family_ids.shape[0] != n_families:
        raise ValueError(f"'family_ids' has {family_ids.shape[0]} along axis 0, but "
            f"'family_centroids' implies n_families == {n_families}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.save_flyer_json_c(
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        expression_vectors,
        family_centroids,
        gene_to_fam,
        is_outlier,
        axis_labels,
        ctypes.byref(ctypes.c_int(axis_labels_strlen)),
        family_ids,
        ctypes.byref(ctypes.c_int(family_ids_strlen)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        gene_species,
        ctypes.byref(ctypes.c_int(gene_species_strlen)),
        gene_types,
        ctypes.byref(ctypes.c_int(gene_types_strlen)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SAVE_FLYER_JSON_ARGUMENTS, _SAVE_FLYER_JSON_ARGUMENT_SOURCES)

    return None
