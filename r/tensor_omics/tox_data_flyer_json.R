# Generated. Do not edit.

#' Write a data set as a flyer JSON file
#'
#' Writes the genes, their expression vectors, their families and the family centroids into
#' a new JSON file, the input of the flyer. Everything is checked before the file is
#' created, and nothing is left behind when a check or a write fails.
#'
#' ### The file
#'
#' One UTF-8 JSON object on a single line, without whitespace, ended by a line feed. With
#' two axes, one family and one gene it reads (broken into lines here):
#'
#' ```
#' {"kind":"tox-flyer","version":1,"tissues":["liver","brain"],
#' "families":[{"family":"F1","gene_indices":[1],"centroid":[5.0000000000000000E-001,-2.0000000000000000E+000]}],
#' "genes":[{"coordinates":[1.0000000000000000E+000,0.0],"id":"g1","family":"F1",
#' "species":"human","is_outlier":false,"type":"ortholog"}]}
#' ```
#'
#' The keys of the object:
#'
#' - `"kind"`: always `"tox-flyer"`;
#' - `"version"`: always `1`;
#' - `"tissues"`: `axis_labels` in order, the name of each coordinate axis;
#' - `"families"`: one object per family, in the order of `family_ids`;
#' - `"genes"`: one object per gene, in input order.
#'
#' The keys of each family:
#'
#' - `"family"`: its id;
#' - `"gene_indices"`: the 1-based positions in `"genes"` of its members, ascending (`[]`
#' for a family without genes);
#' - `"centroid"`: its column of `family_centroids`, one number per axis.
#'
#' The keys of each gene, all six always present:
#'
#' - `"coordinates"`: its column of `expression_vectors`, one number per axis;
#' - `"id"`: its id;
#' - `"family"`: the id of its family, or `null` for a gene in no family (never left out);
#' - `"species"`: its species;
#' - `"is_outlier"`: `true` or `false`;
#' - `"type"`: its type, e.g. ortholog or paralog.
#'
#' Readers ignore keys they do not know, so a key can be added without changing the
#' version. Ids and labels are unique and non-empty. Numbers carry 17 significant digits,
#' which read back as exactly the same double; a zero is written `0.0` or `-0.0`. Strings
#' lose their trailing blanks (so a value that genuinely ends in one cannot be written);
#' `"`, the backslash and the control characters are escaped.
#'
#' ### Empty inputs
#'
#' Zero genes give `"genes":[]` and families whose `"gene_indices"` are all `[]`. Zero
#' families give `"families":[]`, with every gene's `"family"` null. The Python and R
#' bindings cannot pass an empty array of strings yet (issue #221), so there it takes at
#' least one family and one gene. At least one axis is always required.
#'
#' ### Errors
#'
#' Each is reported with the argument it concerns, and the file is not created:
#'
#' - `ERR_EMPTY_INPUT`: no axes;
#' - `ERR_NAN_INF`: a NaN or an infinity in `expression_vectors` or `family_centroids`;
#' - `ERR_INVALID_INPUT`: a `gene_to_fam` entry that is negative or above the number of
#' families; an empty or a repeated entry in `axis_labels`, `family_ids` or `gene_ids`;
#' - `ERR_INVALID_UTF8`: a string that is not valid UTF-8;
#' - `ERR_FILE_OPEN`: `filename` exists already (it is never overwritten) or cannot be
#' created.
#'
#' `ERR_WRITE_DATA` means writing failed, for example on a full disk; the partial file is
#' deleted.
#'
#' Generated from the Fortran procedure \code{tox_data_flyer_json::save_flyer_json}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. Coordinates of each gene, one column per gene
#' @param family_centroids a numeric matrix. Centroid of each family, one column per family
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param is_outlier a logical vector. Whether each gene is an outlier
#' @param axis_labels a character vector. Name of each coordinate axis (e.g. tissues); unique and non-empty
#' @param family_ids a character vector. Id of each family; unique and non-empty
#' @param gene_ids a character vector. Id of each gene; unique and non-empty
#' @param gene_species a character vector. Species of each gene
#' @param gene_types a character vector. Type of each gene, e.g. ortholog or paralog
#' @param filename a string. Path of the JSON file to create; it must not exist yet
#' @return invisibly `NULL`; called for its effect.
#' @export
save_flyer_json <- function(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, gene_ids, gene_species, gene_types, filename) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    is_outlier <- .tox_as_logical_vector(is_outlier, "is_outlier")
    axis_labels <- .tox_as_character(axis_labels, "axis_labels")
    family_ids <- .tox_as_character(family_ids, "family_ids")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    gene_species <- .tox_as_character(gene_species, "gene_species")
    gene_types <- .tox_as_character(gene_types, "gene_types")
    filename <- .tox_as_character(filename, "filename")
    if (dim(family_centroids)[1] != dim(expression_vectors)[1])
        .tox_shape_error("family_centroids", dim(family_centroids)[1], "expression_vectors", dim(expression_vectors)[1])
    if (length(axis_labels) != dim(expression_vectors)[1])
        .tox_shape_error("axis_labels", length(axis_labels), "expression_vectors", dim(expression_vectors)[1])
    if (length(gene_to_fam) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "expression_vectors", dim(expression_vectors)[2])
    if (length(is_outlier) != dim(expression_vectors)[2])
        .tox_shape_error("is_outlier", length(is_outlier), "expression_vectors", dim(expression_vectors)[2])
    if (length(gene_ids) != dim(expression_vectors)[2])
        .tox_shape_error("gene_ids", length(gene_ids), "expression_vectors", dim(expression_vectors)[2])
    if (length(gene_species) != dim(expression_vectors)[2])
        .tox_shape_error("gene_species", length(gene_species), "expression_vectors", dim(expression_vectors)[2])
    if (length(gene_types) != dim(expression_vectors)[2])
        .tox_shape_error("gene_types", length(gene_types), "expression_vectors", dim(expression_vectors)[2])
    if (length(family_ids) != dim(family_centroids)[2])
        .tox_shape_error("family_ids", length(family_ids), "family_centroids", dim(family_centroids)[2])

    .result <- .Call("save_flyer_json_call", expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, gene_ids, gene_species, gene_types, filename)
    .arguments <- c("n_axes", "n_genes", "n_families", "expression_vectors", "family_centroids", "gene_to_fam", "is_outlier", "axis_labels", "family_ids", "gene_ids", "gene_species", "gene_types", "filename", "ierr")
    .sources <- c("expression_vectors", "expression_vectors", "family_centroids", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}
