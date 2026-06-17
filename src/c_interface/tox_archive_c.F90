#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_archive(module)]]
module tox_archive_c
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

    !> summary: C-wrapper for [[tox_archive(module):save_tox_data(subroutine)]]
    !| Save standard tox data
    subroutine save_tox_data_c(zip_filename, zip_filename_strlen, ierr, gene_ids, gene_ids_strlen, n_gene_ids_elements, gene_ids_file, gene_ids_file_strlen, expression, n_expression_elements_dim_1, n_expression_elements_dim_2, expression_file, expression_file_strlen, gene_to_family, n_gene_to_family_elements, gene_to_family_file, gene_to_family_file_strlen, family_ids, family_ids_strlen, n_family_ids_elements, family_ids_file, family_ids_file_strlen, family_centroids, n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2, family_centroids_file, family_centroids_file_strlen, shift_vectors, n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2, shift_vectors_file, shift_vectors_file_strlen) bind(C, name="save_tox_data_c")
        use tox_archive, only: save_tox_data
        use tox_archive
        integer(c_int), intent(in), target :: zip_filename_strlen
            !!  String length of 'zip_filename'
        integer(c_int), intent(in), target :: gene_ids_strlen
            !!  String length of 'gene_ids'
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !!  Size of the 1. dimension/extent of `gene_ids`
        integer(c_int), intent(in), target :: gene_ids_file_strlen
            !!  String length of 'gene_ids_file'
        integer(c_int), intent(in), target :: n_expression_elements_dim_1
            !!  Size of the 1. dimension/extent of `expression`
        integer(c_int), intent(in), target :: n_expression_elements_dim_2
            !!  Size of the 2. dimension/extent of `expression`
        integer(c_int), intent(in), target :: expression_file_strlen
            !!  String length of 'expression_file'
        integer(c_int), intent(in), target :: n_gene_to_family_elements
            !!  Size of the 1. dimension/extent of `gene_to_family`
        integer(c_int), intent(in), target :: gene_to_family_file_strlen
            !!  String length of 'gene_to_family_file'
        integer(c_int), intent(in), target :: family_ids_strlen
            !!  String length of 'family_ids'
        integer(c_int), intent(in), target :: n_family_ids_elements
            !!  Size of the 1. dimension/extent of `family_ids`
        integer(c_int), intent(in), target :: family_ids_file_strlen
            !!  String length of 'family_ids_file'
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !!  Size of the 1. dimension/extent of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !!  Size of the 2. dimension/extent of `family_centroids`
        integer(c_int), intent(in), target :: family_centroids_file_strlen
            !!  String length of 'family_centroids_file'
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_1
            !!  Size of the 1. dimension/extent of `shift_vectors`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_2
            !!  Size of the 2. dimension/extent of `shift_vectors`
        integer(c_int), intent(in), target :: shift_vectors_file_strlen
            !!  String length of 'shift_vectors_file'
        character(len=1, kind=c_char), intent(in), dimension(zip_filename_strlen), target :: zip_filename
            !! Zip filename
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=1, kind=c_char), intent(in), dimension(gene_ids_strlen, n_gene_ids_elements), target :: gene_ids
            !! Gene ids array, will be saved if provided
            !! M_DOC_NO_DEFAULT
        character(len=1, kind=c_char), intent(in), dimension(gene_ids_file_strlen), target :: gene_ids_file
            !! Name of the gene ids file
        real(c_double), intent(in), dimension(n_expression_elements_dim_1, n_expression_elements_dim_2), target :: expression
            !! Expression vectors array, will be saved if provided
        character(len=1, kind=c_char), intent(in), dimension(expression_file_strlen), target :: expression_file
            !! Name of the expression file
        integer(c_int), intent(in), dimension(n_gene_to_family_elements), target :: gene_to_family
            !! Gene to family mapping array, will be saved if provided
        character(len=1, kind=c_char), intent(in), dimension(gene_to_family_file_strlen), target :: gene_to_family_file
            !! Name of the gene to family mapping file
        character(len=1, kind=c_char), intent(in), dimension(family_ids_strlen, n_family_ids_elements), target :: family_ids
            !! Family ids array, will be saved if provided
        character(len=1, kind=c_char), intent(in), dimension(family_ids_file_strlen), target :: family_ids_file
            !! Name of the family ids file
        real(c_double), intent(in), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), target :: family_centroids
            !! Family centroids array, will be saved if provided
        character(len=1, kind=c_char), intent(in), dimension(family_centroids_file_strlen), target :: family_centroids_file
            !! Name of the family centroids file
        real(c_double), intent(in), dimension(n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2), target :: shift_vectors
            !! Shift vectors array, will be saved if provided
        character(len=1, kind=c_char), intent(in), dimension(shift_vectors_file_strlen), target :: shift_vectors_file
            !! Name of the shift vectors file
        character(len=:), allocatable :: zip_filename_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable :: gene_ids_file_f
        character(len=:), allocatable :: expression_file_f
        character(len=:), allocatable :: gene_to_family_file_f
        character(len=:), allocatable, dimension(:) :: family_ids_f
        character(len=:), allocatable :: family_ids_file_f
        character(len=:), allocatable :: family_centroids_file_f
        character(len=:), allocatable :: shift_vectors_file_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(zip_filename)
        M_CHECK_NON_NULL(zip_filename_strlen)
        M_CHECK_NON_NULL(gene_ids)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(gene_ids_file)
        M_CHECK_NON_NULL(gene_ids_file_strlen)
        M_CHECK_NON_NULL(expression)
        M_CHECK_NON_NULL(n_expression_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_elements_dim_2)
        M_CHECK_NON_NULL(expression_file)
        M_CHECK_NON_NULL(expression_file_strlen)
        M_CHECK_NON_NULL(gene_to_family)
        M_CHECK_NON_NULL(n_gene_to_family_elements)
        M_CHECK_NON_NULL(gene_to_family_file)
        M_CHECK_NON_NULL(gene_to_family_file_strlen)
        M_CHECK_NON_NULL(family_ids)
        M_CHECK_NON_NULL(family_ids_strlen)
        M_CHECK_NON_NULL(n_family_ids_elements)
        M_CHECK_NON_NULL(family_ids_file)
        M_CHECK_NON_NULL(family_ids_file_strlen)
        M_CHECK_NON_NULL(family_centroids)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_NON_NULL(family_centroids_file)
        M_CHECK_NON_NULL(family_centroids_file_strlen)
        M_CHECK_NON_NULL(shift_vectors)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_2)
        M_CHECK_NON_NULL(shift_vectors_file)
        M_CHECK_NON_NULL(shift_vectors_file_strlen)
        call c_char_1d_as_string(zip_filename, zip_filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(gene_ids_file, gene_ids_file_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(expression_file, expression_file_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(gene_to_family_file, gene_to_family_file_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(family_ids, family_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(family_ids_file, family_ids_file_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(family_centroids_file, family_centroids_file_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(shift_vectors_file, shift_vectors_file_f, ierr)
        if (is_err(ierr)) return
        call save_tox_data(zip_filename_f, ierr, gene_ids_f, gene_ids_file_f, expression, expression_file_f, gene_to_family, gene_to_family_file_f, family_ids_f, family_ids_file_f, family_centroids, family_centroids_file_f, shift_vectors, shift_vectors_file_f)
    end subroutine save_tox_data_c

end module tox_archive_c
#endif