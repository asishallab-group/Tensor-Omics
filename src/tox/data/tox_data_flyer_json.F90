#include <src/macros.h>

!> AUTHOR_FRANZ_ERIC_SILL
!| Export of a data set to the flyer: a JSON file describing genes, their families and their
!| expression vectors, for viewing them in three dimensions.
!|
!| [[tox_data_flyer_json(module):save_flyer_json(subroutine)]] is the entry point; its
!| documentation is the specification of the file format.
module tox_data_flyer_json
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, set_err_once, is_err, get_err_code, &
                          validate_dimension_size, validate_in_range_int, validate_all_in_range_int, &
                          validate_all_in_range_real, ERR_ALLOC_FAIL, ERR_FILE_OPEN, ERR_INVALID_INPUT, ERR_INVALID_UTF8
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use f42_serde_json_serialize, only: json_writer, json_open, json_close, json_begin_object, json_end_object, &
                                        json_begin_array, json_end_array, json_member, json_null, json_is_valid_utf8
    M_IMPLICIT_NONE

    private
    public :: save_flyer_json

    integer(int32), parameter :: FLYER_FORMAT_VERSION = 1
        !! Value of the file's `"version"` key

contains

    !> M_EXPORT_C
    !| summary: Write a data set as a flyer JSON file
    !| AUTHOR_FRANZ_ERIC_SILL
    !|
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
    !|   for a family without genes);
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
    !|   families; an empty or a repeated entry in `axis_labels`, `family_ids` or `gene_ids`;
    !| - `ERR_INVALID_UTF8`: a string that is not valid UTF-8;
    !| - `ERR_FILE_OPEN`: `filename` exists already (it is never overwritten) or cannot be
    !|   created.
    !|
    !| `ERR_WRITE_DATA` means writing failed, for example on a full disk; the partial file is
    !| deleted.
    subroutine save_flyer_json(n_axes, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                               is_outlier, axis_labels, family_ids, gene_ids, gene_species, gene_types, filename, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of coordinate axes (e.g. tissues)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! Coordinates of each gene, one column per gene
        real(real64), dimension(n_axes, n_families), intent(in) :: family_centroids
            !! Centroid of each family, one column per family
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        logical(c_bool), dimension(n_genes), intent(in) :: is_outlier
            !! Whether each gene is an outlier
        character(len=*), dimension(n_axes), intent(in) :: axis_labels
            !! Name of each coordinate axis (e.g. tissues); unique and non-empty
        character(len=*), dimension(n_families), intent(in) :: family_ids
            !! Id of each family; unique and non-empty
        character(len=*), dimension(n_genes), intent(in) :: gene_ids
            !! Id of each gene; unique and non-empty
        character(len=*), dimension(n_genes), intent(in) :: gene_species
            !! Species of each gene
        character(len=*), dimension(n_genes), intent(in) :: gene_types
            !! Type of each gene, e.g. ortholog or paralog
        character(len=*), intent(in) :: filename
            !! Path of the JSON file to create; it must not exist yet
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), allocatable :: sort_perm(:), family_start(:), family_members(:)
        integer(int32) :: i_gene, i_family
        type(json_writer) :: writer

        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=1_int32)
        call validate_in_range_int(n_genes, ierr, min=0_int32, arg_pos=2_int32)
        call validate_in_range_int(n_families, ierr, min=0_int32, arg_pos=3_int32)
        if (is_err(ierr)) return

        ! One permutation serves every duplicate check, so it is sized for the longest list.
        M_ALLOCATE(sort_perm(max(n_axes, n_genes, n_families)))

        ! column by column, so no element count has to fit an int32
        do i_gene = 1, n_genes
            call validate_all_in_range_real(expression_vectors(:, i_gene), n_axes, ierr, arg_pos=4_int32)
        end do
        do i_family = 1, n_families
            call validate_all_in_range_real(family_centroids(:, i_family), n_axes, ierr, arg_pos=5_int32)
        end do
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, min=1_int32, max=n_families, &
                                       sentinel=M_GENE_TO_FAM_SENTINEL, arg_pos=6_int32)
        call validate_labels(axis_labels, .true._c_bool, sort_perm, 8_int32, ierr)
        call validate_labels(family_ids, .true._c_bool, sort_perm, 9_int32, ierr)
        call validate_labels(gene_ids, .true._c_bool, sort_perm, 10_int32, ierr)
        call validate_labels(gene_species, .false._c_bool, sort_perm, 11_int32, ierr)
        call validate_labels(gene_types, .false._c_bool, sort_perm, 12_int32, ierr)
        if (is_err(ierr)) return

        call group_genes_by_family(gene_to_fam, n_families, family_start, family_members, ierr)
        if (is_err(ierr)) return

        call json_open(writer, filename, ierr)
        if (is_err(ierr)) then
            if (get_err_code(ierr) == ERR_FILE_OPEN) call set_err(ierr, ERR_FILE_OPEN, 13_int32)
            return
        end if

        call json_begin_object(writer)
        call json_member(writer, 'kind', 'tox-flyer')
        call json_member(writer, 'version', FLYER_FORMAT_VERSION)
        call json_member(writer, 'tissues', axis_labels)

        call json_begin_array(writer, 'families')
        do i_family = 1, n_families
            call json_begin_object(writer)
            call json_member(writer, 'family', family_ids(i_family))
            call json_member(writer, 'gene_indices', &
                             family_members(family_start(i_family):family_start(i_family + 1) - 1))
            call json_member(writer, 'centroid', family_centroids(:, i_family))
            call json_end_object(writer)
        end do
        call json_end_array(writer)

        call json_begin_array(writer, 'genes')
        do i_gene = 1, n_genes
            call json_begin_object(writer)
            call json_member(writer, 'coordinates', expression_vectors(:, i_gene))
            call json_member(writer, 'id', gene_ids(i_gene))
            if (gene_to_fam(i_gene) /= M_GENE_TO_FAM_SENTINEL) then
                call json_member(writer, 'family', family_ids(gene_to_fam(i_gene)))
            else
                call json_null(writer, 'family')
            end if
            call json_member(writer, 'species', gene_species(i_gene))
            call json_member(writer, 'is_outlier', is_outlier(i_gene))
            call json_member(writer, 'type', gene_types(i_gene))
            call json_end_object(writer)
        end do
        call json_end_array(writer)
        call json_end_object(writer)

        ! The writer fails only on writing, everything else was checked above; its error is
        ! about the file as a whole, not about one argument.
        call json_close(writer, ierr)
    end subroutine save_flyer_json

    !> The genes of each family, in ascending order, as compressed rows: family `f` owns
    !| `family_members(family_start(f):family_start(f + 1) - 1)`.
    pure subroutine group_genes_by_family(gene_to_fam, n_families, family_start, family_members, ierr)
        integer(int32), intent(in) :: gene_to_fam(:)
            !! Family of each gene, already validated; `M_GENE_TO_FAM_SENTINEL` for none
        integer(int32), intent(in) :: n_families
            !! Number of families
        integer(int32), allocatable, intent(out) :: family_start(:)
            !! Where each family's genes start in `family_members`, and one past the last family's end
        integer(int32), allocatable, intent(out) :: family_members(:)
            !! The genes of every family, family after family
        integer(int32), intent(out) :: ierr
            !! Error code: `ERR_OK` or `ERR_ALLOC_FAIL`

        integer(int32) :: i_gene, i_family, family_of_gene

        call set_ok(ierr)
        M_ALLOCATE(family_start(n_families + 1))
        family_start = 0
        do i_gene = 1, size(gene_to_fam, kind=int32)
            family_of_gene = gene_to_fam(i_gene)
            if (family_of_gene /= M_GENE_TO_FAM_SENTINEL) &
                family_start(family_of_gene + 1) = family_start(family_of_gene + 1) + 1
        end do
        family_start(1) = 1
        do i_family = 2, n_families + 1
            family_start(i_family) = family_start(i_family) + family_start(i_family - 1)
        end do
        M_ALLOCATE(family_members(family_start(n_families + 1) - 1))
        ! Filling advances family_start(f) to the start of family f + 1 ...
        do i_gene = 1, size(gene_to_fam, kind=int32)
            family_of_gene = gene_to_fam(i_gene)
            if (family_of_gene /= M_GENE_TO_FAM_SENTINEL) then
                family_members(family_start(family_of_gene)) = i_gene
                family_start(family_of_gene) = family_start(family_of_gene) + 1
            end if
        end do
        ! ... so every start moves back by one family, walking down to leave the unread ones intact
        ! (family_start(n_families + 1), one past the end, is already right).
        do i_family = n_families, 2, -1
            family_start(i_family) = family_start(i_family - 1)
        end do
        family_start(1) = 1
    end subroutine group_genes_by_family

    !> Check a list of strings for the file: valid UTF-8, and for identifiers also non-empty
    !| and unique.
    !|
    !| Duplicates are found by sorting a permutation and comparing neighbours, which compares
    !| the strings without their trailing blanks, as they are written.
    subroutine validate_labels(labels, are_identifiers, sort_perm, arg_pos, ierr)
        character(len=*), intent(in), contiguous :: labels(:)
            !! The strings
        logical(c_bool), intent(in) :: are_identifiers
            !! Whether they must also be non-empty and unique
        integer(int32), intent(inout), contiguous :: sort_perm(:)
            !! Work array of at least `size(labels)` elements
        integer(int32), intent(in) :: arg_pos
            !! Position of the argument holding them
        integer(int32), intent(inout) :: ierr
            !! Error code; only a first error is recorded

        integer(int32) :: n_labels, i_label

        n_labels = size(labels, kind=int32)
        do i_label = 1, n_labels
            if (are_identifiers .and. len_trim(labels(i_label)) == 0) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
                return
            end if
            if (.not. json_is_valid_utf8(labels(i_label))) then
                call set_err_once(ierr, ERR_INVALID_UTF8, arg_pos)
                return
            end if
        end do
        if (.not. are_identifiers .or. n_labels < 2) return

        call init_perm(sort_perm(1:n_labels))
        call sort_array_heapsort(labels, sort_perm(1:n_labels))
        do i_label = 2, n_labels
            if (labels(sort_perm(i_label)) == labels(sort_perm(i_label - 1))) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
                return
            end if
        end do
    end subroutine validate_labels

end module tox_data_flyer_json
