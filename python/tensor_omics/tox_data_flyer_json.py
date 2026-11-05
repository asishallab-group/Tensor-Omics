r"""tox_data_flyer_json

Export of a data set to the flyer: a JSON file describing genes, their families and their
expression vectors, for viewing them in three dimensions.

:func:`tensor_omics.save_flyer_json` is the entry point; its
documentation is the specification of the file format.

Python binding, generated from tox_data_flyer_json. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.save_flyer_json_c.restype = None
_lib.save_flyer_json_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    family_ids : sequence of str, of length n_families
        Id of each family; unique and non-empty
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_ids : sequence of str, of length n_genes
        Id of each gene; unique and non-empty
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_species : sequence of str, of length n_genes
        Species of each gene
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    gene_types : sequence of str, of length n_genes
        Type of each gene, e.g. ortholog or paralog
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    filename : str
        Path of the JSON file to create; it must not exist yet
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

    Returns
    -------
    None

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_flyer_json::save_flyer_json`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    family_centroids = _contract.reals(family_centroids, "family_centroids", ndim=2, order="F")
    gene_to_fam = _contract.integers(gene_to_fam, "gene_to_fam", ndim=1)
    is_outlier = _contract.logicals(is_outlier, "is_outlier", ndim=1)
    _axis_labels_bytes = _contract.strings(axis_labels, "axis_labels")
    _axis_labels_width = max(map(len, _axis_labels_bytes), default=0) or 1
    axis_labels = _np.asarray([_b.ljust(_axis_labels_width) for _b in _axis_labels_bytes], dtype="S")
    _family_ids_bytes = _contract.strings(family_ids, "family_ids")
    _family_ids_width = max(map(len, _family_ids_bytes), default=0) or 1
    family_ids = _np.asarray([_b.ljust(_family_ids_width) for _b in _family_ids_bytes], dtype="S")
    _gene_ids_bytes = _contract.strings(gene_ids, "gene_ids")
    _gene_ids_width = max(map(len, _gene_ids_bytes), default=0) or 1
    gene_ids = _np.asarray([_b.ljust(_gene_ids_width) for _b in _gene_ids_bytes], dtype="S")
    _gene_species_bytes = _contract.strings(gene_species, "gene_species")
    _gene_species_width = max(map(len, _gene_species_bytes), default=0) or 1
    gene_species = _np.asarray([_b.ljust(_gene_species_width) for _b in _gene_species_bytes], dtype="S")
    _gene_types_bytes = _contract.strings(gene_types, "gene_types")
    _gene_types_width = max(map(len, _gene_types_bytes), default=0) or 1
    gene_types = _np.asarray([_b.ljust(_gene_types_width) for _b in _gene_types_bytes], dtype="S")
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")

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

    # Fortran cannot check that an array has the shape it is declared with; this can
    if family_centroids.shape[0] != n_axes:
        raise _ToxArgumentError(f"'family_centroids' has {family_centroids.shape[0]} along axis 0, but 'expression_vectors' implies n_axes == {n_axes}; 'family_centroids' is declared with shape (n_axes, n_families)", argument="family_centroids")
    if gene_to_fam.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_to_fam' is declared with shape (n_genes)", argument="gene_to_fam")
    if is_outlier.shape[0] != n_genes:
        raise _ToxArgumentError(f"'is_outlier' has {is_outlier.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'is_outlier' is declared with shape (n_genes)", argument="is_outlier")
    if axis_labels.shape[0] != n_axes:
        raise _ToxArgumentError(f"'axis_labels' has {axis_labels.shape[0]} along axis 0, but 'expression_vectors' implies n_axes == {n_axes}; 'axis_labels' is declared with shape (n_axes)", argument="axis_labels")
    if family_ids.shape[0] != n_families:
        raise _ToxArgumentError(f"'family_ids' has {family_ids.shape[0]} along axis 0, but 'family_centroids' implies n_families == {n_families}; 'family_ids' is declared with shape (n_families)", argument="family_ids")
    if gene_ids.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_ids' has {gene_ids.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_ids' is declared with shape (n_genes)", argument="gene_ids")
    if gene_species.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_species' has {gene_species.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_species' is declared with shape (n_genes)", argument="gene_species")
    if gene_types.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_types' has {gene_types.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_types' is declared with shape (n_genes)", argument="gene_types")

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.save_flyer_json_c(
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_families)),
        expression_vectors,
        family_centroids,
        gene_to_fam,
        is_outlier,
        axis_labels,
        _ctypes.byref(_ctypes.c_int(axis_labels_strlen)),
        family_ids,
        _ctypes.byref(_ctypes.c_int(family_ids_strlen)),
        gene_ids,
        _ctypes.byref(_ctypes.c_int(gene_ids_strlen)),
        gene_species,
        _ctypes.byref(_ctypes.c_int(gene_species_strlen)),
        gene_types,
        _ctypes.byref(_ctypes.c_int(gene_types_strlen)),
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _SAVE_FLYER_JSON_ARGUMENTS, _SAVE_FLYER_JSON_ARGUMENT_SOURCES)

    return None
