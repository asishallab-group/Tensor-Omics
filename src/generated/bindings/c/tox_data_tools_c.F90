#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_tools(module)]]
!| Parsers for the plain-text input formats TensorOmics data sets are built from (gene-expression
!| TSV/CSV files, OrthoFinder-style family files), plus small array-filtering helpers.
!|
!| These populate the raw `gene_ids` / `expression` / `gene_to_family` / `family_ids` arrays that
!| [[tox_data_archive(module):save_tox_data(subroutine)]] later persists; unlike the archive/serde
!| layers, everything here works from delimited text rather than the library's binary array format.
module tox_data_tools_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string, string_as_c_char_2d
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: read_expression_vectors_tsv_c
    public :: read_gene_ids_from_tsv_file_c
    public :: read_orthofinder_file_c
    public :: get_unassigned_mask_c

contains

    !> summary: C-wrapper for [[tox_data_tools(module):read_expression_vectors_tsv(subroutine)]]
    subroutine read_expression_vectors_tsv_c(&
            file_list,&
            file_list_strlen,&
            n_file_list_elements,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            expression_vectors,&
            n_expression_vectors_elements_dim_1,&
            n_expression_vectors_elements_dim_2,&
            n_header_rows,&
            gene_col,&
            value_cols,&
            n_value_cols_elements,&
            start_row,&
            ierr,&
            delimiter&
        ) bind(C, name="read_expression_vectors_tsv_c")
        use tox_data_tools, only: read_expression_vectors_tsv

        integer(c_int), intent(in), target :: file_list_strlen
            !! length of the strings in `file_list`
        integer(c_int), intent(in), target :: n_file_list_elements
            !! number of elements in `file_list`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !! 1. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !! 2. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_value_cols_elements
            !! number of elements in `value_cols`
        character(len=1, kind=c_char), dimension(file_list_strlen, n_file_list_elements), intent(in), target :: file_list
            !! List of files to read from
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(in), target :: gene_ids
            !! Array of gene IDS
        real(c_double), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), intent(inout), target :: expression_vectors
            !! Array of expression vectors
        integer(c_int), intent(in), target :: n_header_rows
            !! Number of header rows to skip
        integer(c_int), intent(in), target :: gene_col
            !! Index of column with gene_ids
        integer(c_int), dimension(n_value_cols_elements), intent(in), target :: value_cols
            !! Indicies of columns containing values
        integer(c_int), intent(in), target :: start_row
            !! Row in the expression vectors to start in
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=1, kind=c_char), dimension(1), intent(in), target :: delimiter
            !! optional delimiter
            !! The default value is `char(9)`.
        character(len=:), allocatable, dimension(:) :: file_list_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=1), allocatable :: delimiter_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(file_list_strlen)
        M_CHECK_NON_NULL(n_file_list_elements)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_header_rows)
        M_CHECK_NON_NULL(gene_col)
        M_CHECK_NON_NULL(n_value_cols_elements)
        M_CHECK_NON_NULL(start_row)
        M_CHECK_ARRAY_NON_NULL(file_list, file_list_strlen * n_file_list_elements)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_gene_ids_elements)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_expression_vectors_elements_dim_1 * n_expression_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(value_cols, n_value_cols_elements)
        M_CHECK_ARRAY_NON_NULL(delimiter, 1)

        call c_char_2d_as_string(file_list, file_list_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        M_ALLOCATE(delimiter_f)
        block
            character(len=:), allocatable :: converted
            call c_char_1d_as_string(delimiter, converted, ierr)
            if (is_err(ierr)) return
            delimiter_f = converted
        end block

        call read_expression_vectors_tsv(&
            file_list = file_list_f,&
            gene_ids = gene_ids_f,&
            expression_vectors = expression_vectors,&
            n_header_rows = n_header_rows,&
            gene_col = gene_col,&
            value_cols = value_cols,&
            start_row = start_row,&
            ierr = ierr,&
            delimiter = delimiter_f&
        )
    end subroutine read_expression_vectors_tsv_c

    !> summary: C-wrapper for [[tox_data_tools(module):read_gene_ids_from_tsv_file(subroutine)]]
    subroutine read_gene_ids_from_tsv_file_c(&
            filename,&
            filename_strlen,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            n_header_rows,&
            gene_col,&
            ierr&
        ) bind(C, name="read_gene_ids_from_tsv_file_c")
        use tox_data_tools, only: read_gene_ids_from_tsv_file

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(out), target :: gene_ids
            !! gene ids array
        integer(c_int), intent(in), target :: n_header_rows
            !! number of headers to skip
        integer(c_int), intent(in), target :: gene_col
            !! Index of the column containing gene ids
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: filename_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(n_header_rows)
        M_CHECK_NON_NULL(gene_col)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_gene_ids_elements)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        M_ALLOCATE(character(len=gene_ids_strlen) :: gene_ids_f(n_gene_ids_elements))

        call read_gene_ids_from_tsv_file(&
            filename = filename_f,&
            gene_ids = gene_ids_f,&
            n_header_rows = n_header_rows,&
            gene_col = gene_col,&
            ierr = ierr&
        )

        call string_as_c_char_2d(gene_ids_f, gene_ids)
    end subroutine read_gene_ids_from_tsv_file_c

    !> summary: C-wrapper for [[tox_data_tools(module):read_orthofinder_file(subroutine)]]
    subroutine read_orthofinder_file_c(&
            filename,&
            filename_strlen,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            family_ids,&
            family_ids_strlen,&
            n_family_ids_elements,&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            ierr&
        ) bind(C, name="read_orthofinder_file_c")
        use tox_data_tools, only: read_orthofinder_file

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        integer(c_int), intent(in), target :: family_ids_strlen
            !! length of the strings in `family_ids`
        integer(c_int), intent(in), target :: n_family_ids_elements
            !! number of elements in `family_ids`
        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(in), target :: gene_ids
            !! gene ids array
        character(len=1, kind=c_char), dimension(family_ids_strlen, n_family_ids_elements), intent(out), target :: family_ids
            !! family ids array
        integer(c_int), dimension(n_gene_to_fam_elements), intent(out), target :: gene_to_fam
            !! gene to family mapping
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: filename_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable, dimension(:) :: family_ids_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(family_ids_strlen)
        M_CHECK_NON_NULL(n_family_ids_elements)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_gene_ids_elements)
        M_CHECK_ARRAY_NON_NULL(family_ids, family_ids_strlen * n_family_ids_elements)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        M_ALLOCATE(character(len=family_ids_strlen) :: family_ids_f(n_family_ids_elements))

        call read_orthofinder_file(&
            filename = filename_f,&
            gene_ids = gene_ids_f,&
            family_ids = family_ids_f,&
            gene_to_fam = gene_to_fam,&
            ierr = ierr&
        )

        call string_as_c_char_2d(family_ids_f, family_ids)
    end subroutine read_orthofinder_file_c

    !> summary: C-wrapper for [[tox_data_tools(module):get_unassigned_mask(subroutine)]]
    subroutine get_unassigned_mask_c(&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            mask,&
            n_genes_kept,&
            ierr&
        ) bind(C, name="get_unassigned_mask_c")
        use tox_data_tools, only: get_unassigned_mask

        integer(c_int), dimension(n_gene_to_fam_elements), intent(in), target :: gene_to_fam
            !! gene to family mapping
        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        logical(c_bool), dimension(size(gene_to_fam)), intent(out), target :: mask
            !! mask for mapping
        integer(c_int), intent(out), target :: n_genes_kept
            !! number of genes kept
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(:), allocatable :: mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_NON_NULL(n_genes_kept)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)
        M_CHECK_ARRAY_NON_NULL(mask, (size(gene_to_fam)))

        M_ALLOCATE(mask_f(size(gene_to_fam)))

        call get_unassigned_mask(&
            gene_to_fam = gene_to_fam,&
            mask = mask_f,&
            n_genes_kept = n_genes_kept&
        )

        mask = mask_f
    end subroutine get_unassigned_mask_c

end module tox_data_tools_c
#endif
