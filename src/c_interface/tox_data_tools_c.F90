#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_data_tools(module)]]
module tox_data_tools_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[tox_data_tools(module):read_expression_vectors_tsv(subroutine)]]
    !| Read expression vectors from csv/tsv files
    subroutine read_expression_vectors_tsv_c(file_list, file_list_strlen, n_file_list_elements, gene_ids, gene_ids_strlen, n_gene_ids_elements, expression_vectors, n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2, n_header_rows, gene_col, value_cols, n_value_cols_elements, start_row, ierr, delimiter) bind(C, name="read_expression_vectors_tsv_c")
        use tox_data_tools, only: read_expression_vectors_tsv
        use tox_data_tools
        integer(c_int), intent(in), target :: file_list_strlen
            !!  String length of 'file_list'
        integer(c_int), intent(in), target :: n_file_list_elements
            !!  Size of the 1. dimension/extent of `file_list`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !!  String length of 'gene_ids'
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !!  Size of the 1. dimension/extent of `gene_ids`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !!  Size of the 1. dimension/extent of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !!  Size of the 2. dimension/extent of `expression_vectors`
        integer(c_int), intent(in), target :: n_value_cols_elements
            !!  Size of the 1. dimension/extent of `value_cols`
        character(len=1, kind=c_char), intent(in), dimension(file_list_strlen, n_file_list_elements), target :: file_list
            !! List of files to read from
        character(len=1, kind=c_char), intent(in), dimension(gene_ids_strlen, n_gene_ids_elements), target :: gene_ids
            !! Array of gene IDS
        real(c_double), intent(inout), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), target :: expression_vectors
            !! Array of expression vectors
        integer(c_int), intent(in), target :: n_header_rows
            !! Number of header rows to skip
        integer(c_int), intent(in), target :: gene_col
            !! Index of column with gene_ids
        integer(c_int), intent(in), dimension(n_value_cols_elements), target :: value_cols
            !! Indicies of columns containing values
        integer(c_int), intent(in), target :: start_row
            !! Row in the expression vectors to start in
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=1, kind=c_char), intent(in), dimension(1), target :: delimiter
            !! optional delimiter
            !! M_DOC_DEFAULT('\t')
        character(len=:), allocatable, dimension(:) :: file_list_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable :: delimiter_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(file_list)
        M_CHECK_NON_NULL(file_list_strlen)
        M_CHECK_NON_NULL(n_file_list_elements)
        M_CHECK_NON_NULL(gene_ids)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(expression_vectors)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_header_rows)
        M_CHECK_NON_NULL(gene_col)
        M_CHECK_NON_NULL(value_cols)
        M_CHECK_NON_NULL(n_value_cols_elements)
        M_CHECK_NON_NULL(start_row)
        M_CHECK_NON_NULL(delimiter)
        call c_char_2d_as_string(file_list, file_list_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(delimiter, delimiter_f, ierr)
        if (is_err(ierr)) return
        call read_expression_vectors_tsv(file_list_f, gene_ids_f, expression_vectors, n_header_rows, gene_col, value_cols, start_row, ierr, delimiter_f)
    end subroutine read_expression_vectors_tsv_c

end module tox_data_tools_c
#endif