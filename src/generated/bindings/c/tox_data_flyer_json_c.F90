#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_flyer_json(module)]]
!| Export of a data set to the flyer: a JSON file describing genes, their families and their
!| expression vectors, for viewing them in three dimensions.
!|
!| [[tox_data_flyer_json(module):save_flyer_json(subroutine)]] is the entry point; its
!| documentation is the specification of the file format.
module tox_data_flyer_json_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_f_pointer, c_int
    use, intrinsic :: iso_c_binding, only: c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_EMPTY_INPUT
    M_IMPLICIT_NONE
    private

    public :: save_flyer_json_c

contains

    !> summary: C-wrapper for [[tox_data_flyer_json(module):save_flyer_json(subroutine)]]
    !| Writes the genes, their expression vectors, their families and the family centroids into
    !| a new JSON file, the input of the flyer. Everything is checked before the file is
    !| created, and nothing is left behind when a check or a write fails.
    !|
    !| ### The file
    !|
    !| One UTF-8 JSON object on a single line, without whitespace, ended by a line feed. With
    !| two axes, one family and one gene it reads (broken into lines here):
    !|
    !| ```
    !| {"kind":"tox-flyer","version":1,"tissues":["liver","brain"],
    !| "families":[{"family":"F1","gene_indices":[1],"centroid":[5.0000000000000000E-001,-2.0000000000000000E+000]}],
    !| "genes":[{"coordinates":[1.0000000000000000E+000,0.0],"id":"g1","family":"F1",
    !| "species":"human","is_outlier":false,"type":"ortholog"}]}
    !| ```
    !|
    !| The keys of the object:
    !|
    !| - `"kind"`: always `"tox-flyer"`;
    !| - `"version"`: always `1`;
    !| - `"tissues"`: `axis_labels` in order, the name of each coordinate axis;
    !| - `"families"`: one object per family, in the order of `family_ids`;
    !| - `"genes"`: one object per gene, in input order.
    !|
    !| The keys of each family:
    !|
    !| - `"family"`: its id;
    !| - `"gene_indices"`: the 1-based positions in `"genes"` of its members, ascending (`[]`
    !| for a family without genes);
    !| - `"centroid"`: its column of `family_centroids`, one number per axis.
    !|
    !| The keys of each gene, all six always present:
    !|
    !| - `"coordinates"`: its column of `expression_vectors`, one number per axis;
    !| - `"id"`: its id;
    !| - `"family"`: the id of its family, or `null` for a gene in no family (never left out);
    !| - `"species"`: its species;
    !| - `"is_outlier"`: `true` or `false`;
    !| - `"type"`: its type, e.g. ortholog or paralog.
    !|
    !| Readers ignore keys they do not know, so a key can be added without changing the
    !| version. Ids and labels are unique and non-empty. Numbers carry 17 significant digits,
    !| which read back as exactly the same double; a zero is written `0.0` or `-0.0`. Strings
    !| lose their trailing blanks (so a value that genuinely ends in one cannot be written);
    !| `"`, the backslash and the control characters are escaped.
    !|
    !| ### Empty inputs
    !|
    !| Zero genes give `"genes":[]` and families whose `"gene_indices"` are all `[]`. Zero
    !| families give `"families":[]`, with every gene's `"family"` null. The Python and R
    !| bindings cannot pass an empty array of strings yet (issue #221), so there it takes at
    !| least one family and one gene. At least one axis is always required.
    !|
    !| ### Errors
    !|
    !| Each is reported with the argument it concerns, and the file is not created:
    !|
    !| - `ERR_EMPTY_INPUT`: no axes;
    !| - `ERR_NAN_INF`: a NaN or an infinity in `expression_vectors` or `family_centroids`;
    !| - `ERR_INVALID_INPUT`: a `gene_to_fam` entry that is negative or above the number of
    !| families; an empty or a repeated entry in `axis_labels`, `family_ids` or `gene_ids`;
    !| - `ERR_INVALID_UTF8`: a string that is not valid UTF-8;
    !| - `ERR_FILE_OPEN`: `filename` exists already (it is never overwritten) or cannot be
    !| created.
    !|
    !| `ERR_WRITE_DATA` means writing failed, for example on a full disk; the partial file is
    !| deleted.
    subroutine save_flyer_json_c(&
            n_axes,&
            n_genes,&
            n_families,&
            expression_vectors,&
            family_centroids,&
            gene_to_fam,&
            is_outlier,&
            axis_labels,&
            axis_labels_strlen,&
            family_ids,&
            family_ids_strlen,&
            gene_ids,&
            gene_ids_strlen,&
            gene_species,&
            gene_species_strlen,&
            gene_types,&
            gene_types_strlen,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="save_flyer_json_c")
        use tox_data_flyer_json, only: save_flyer_json

        integer(c_int), intent(in), target :: n_axes
            !! Number of coordinate axes (e.g. tissues)
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        integer(c_int), intent(in), target :: axis_labels_strlen
            !! length of the strings in `axis_labels`
        integer(c_int), intent(in), target :: family_ids_strlen
            !! length of the strings in `family_ids`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: gene_species_strlen
            !! length of the strings in `gene_species`
        integer(c_int), intent(in), target :: gene_types_strlen
            !! length of the strings in `gene_types`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! Coordinates of each gene, one column per gene
        real(c_double), dimension(n_axes, n_families), intent(in), target :: family_centroids
            !! Centroid of each family, one column per family
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        logical(c_bool), dimension(n_genes), intent(in), target :: is_outlier
            !! Whether each gene is an outlier
        character(len=1, kind=c_char), dimension(axis_labels_strlen, n_axes), intent(in), target :: axis_labels
            !! Name of each coordinate axis (e.g. tissues); unique and non-empty
        character(len=1, kind=c_char), dimension(family_ids_strlen, n_families), intent(in), target :: family_ids
            !! Id of each family; unique and non-empty
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_genes), intent(in), target :: gene_ids
            !! Id of each gene; unique and non-empty
        character(len=1, kind=c_char), dimension(gene_species_strlen, n_genes), intent(in), target :: gene_species
            !! Species of each gene
        character(len=1, kind=c_char), dimension(gene_types_strlen, n_genes), intent(in), target :: gene_types
            !! Type of each gene, e.g. ortholog or paralog
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Path of the JSON file to create; it must not exist yet
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=axis_labels_strlen), pointer, dimension(:) :: axis_labels_f
        character(len=family_ids_strlen), pointer, dimension(:) :: family_ids_f
        character(len=gene_ids_strlen), pointer, dimension(:) :: gene_ids_f
        character(len=gene_species_strlen), pointer, dimension(:) :: gene_species_f
        character(len=gene_types_strlen), pointer, dimension(:) :: gene_types_f
        character(len=filename_strlen), pointer :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(axis_labels_strlen)
        M_CHECK_NON_NULL(family_ids_strlen)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(gene_species_strlen)
        M_CHECK_NON_NULL(gene_types_strlen)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(axis_labels, axis_labels_strlen * n_axes)
        M_CHECK_ARRAY_NON_NULL(family_ids, family_ids_strlen * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_species, gene_species_strlen * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_types, gene_types_strlen * n_genes)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        M_CHECK_CHARACTER_VIEW(axis_labels, axis_labels_strlen * n_axes)
        call c_f_pointer(c_loc(axis_labels), axis_labels_f, [n_axes])
        M_CHECK_CHARACTER_VIEW(family_ids, family_ids_strlen * n_families)
        call c_f_pointer(c_loc(family_ids), family_ids_f, [n_families])
        M_CHECK_CHARACTER_VIEW(gene_ids, gene_ids_strlen * n_genes)
        call c_f_pointer(c_loc(gene_ids), gene_ids_f, [n_genes])
        M_CHECK_CHARACTER_VIEW(gene_species, gene_species_strlen * n_genes)
        call c_f_pointer(c_loc(gene_species), gene_species_f, [n_genes])
        M_CHECK_CHARACTER_VIEW(gene_types, gene_types_strlen * n_genes)
        call c_f_pointer(c_loc(gene_types), gene_types_f, [n_genes])
        M_CHECK_CHARACTER_VIEW(filename, filename_strlen)
        call c_f_pointer(c_loc(filename), filename_f)

        call save_flyer_json(&
            n_axes = n_axes,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            gene_to_fam = gene_to_fam,&
            is_outlier = is_outlier,&
            axis_labels = axis_labels_f,&
            family_ids = family_ids_f,&
            gene_ids = gene_ids_f,&
            gene_species = gene_species_f,&
            gene_types = gene_types_f,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine save_flyer_json_c

end module tox_data_flyer_json_c
#endif
