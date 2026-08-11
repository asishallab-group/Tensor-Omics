# Generated. Do not edit.

#' Serializes tox related data to JSON, compatible with the tox_flyer
#'
#' Generated from the Fortran procedure \code{tox_flyer_json::serialize_tox_data_as_flyer_json}, whose argument names
#' are the ones an error message reports.
#'
#' @param filename a string. Name of the file to write the output to
#' @param tissues a character vector. Tissue identifiers
#' @param family_ids a character vector. Family identifiers
#' @param centroids a numeric matrix. Centroid data
#' @param gene_ids a character vector. Gene identifiers
#' @param genes a numeric matrix. Gene data
#' @param gene_to_fam a integer vector. Gene index to Family index mapping
#' @param sorted_gene_to_fam_perm a integer vector. Permutation vector that sorts `gene_to_fam`
#' @param gene_outliers a logical vector. Specifies if a gene is an outlier
#' @param gene_species a character vector. Species name per gene
#' @param gene_types a character vector. Gene type string (ortholog/paralog)
#' @return invisibly `NULL`; called for its effect.
#' @export
serialize_tox_data_as_flyer_json <- function(filename, tissues, family_ids, centroids, gene_ids, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types) {
    filename <- .tox_as_character(filename, "filename")
    tissues <- .tox_as_character(tissues, "tissues")
    family_ids <- .tox_as_character(family_ids, "family_ids")
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    genes <- .tox_as_double_matrix(genes, "genes")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    sorted_gene_to_fam_perm <- .tox_as_integer_vector(sorted_gene_to_fam_perm, "sorted_gene_to_fam_perm")
    gene_outliers <- .tox_as_logical(gene_outliers, "gene_outliers")
    gene_species <- .tox_as_character(gene_species, "gene_species")
    gene_types <- .tox_as_character(gene_types, "gene_types")
    if (dim(centroids)[1] != length(tissues))
        .tox_shape_error("centroids", dim(centroids)[1], "tissues", length(tissues))
    if (dim(genes)[1] != length(tissues))
        .tox_shape_error("genes", dim(genes)[1], "tissues", length(tissues))
    if (dim(centroids)[2] != length(family_ids))
        .tox_shape_error("centroids", dim(centroids)[2], "family_ids", length(family_ids))
    if (dim(genes)[2] != length(gene_ids))
        .tox_shape_error("genes", dim(genes)[2], "gene_ids", length(gene_ids))
    if (length(gene_to_fam) != length(gene_ids))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "gene_ids", length(gene_ids))
    if (length(sorted_gene_to_fam_perm) != length(gene_ids))
        .tox_shape_error("sorted_gene_to_fam_perm", length(sorted_gene_to_fam_perm), "gene_ids", length(gene_ids))
    if (length(gene_outliers) != length(gene_ids))
        .tox_shape_error("gene_outliers", length(gene_outliers), "gene_ids", length(gene_ids))
    if (length(gene_species) != length(gene_ids))
        .tox_shape_error("gene_species", length(gene_species), "gene_ids", length(gene_ids))
    if (length(gene_types) != length(gene_ids))
        .tox_shape_error("gene_types", length(gene_types), "gene_ids", length(gene_ids))

    .result <- .Call("serialize_tox_data_as_flyer_json_call", filename, tissues, family_ids, centroids, gene_ids, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types)
    .arguments <- c("filename", "tissues", "n_tissues", "family_ids", "n_families", "centroids", "gene_ids", "n_genes", "genes", "gene_to_fam", "sorted_gene_to_fam_perm", "gene_outliers", "gene_species", "gene_types", "ierr")
    .sources <- c(NA_character_, NA_character_, "tissues", NA_character_, "family_ids", NA_character_, NA_character_, "gene_ids", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    invisible(NULL)
}
