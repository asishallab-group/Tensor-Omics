#include <src/macros.h>

!> Module for computing expression centroids of gene families.
!|
!| This module contains the core scientific kernel. The C and R interface
!| wrappers are defined outside the module for compatibility.
module tox_gene_centroids
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, validate_dimension_size, validate_in_range_int, validate_all_in_range_int
    implicit none

#define CM_MODE_GROUP_ORTHOLOGS 0_int32
#define CM_MODE_GROUP_ALL 1_int32

    integer(int32), parameter, public :: MODE_GROUP_ORTHOLOGS = CM_MODE_GROUP_ORTHOLOGS
        !! Mode code for grouping by orthologs in [[tox_gene_centroids(module):group_centroid(subroutine)]]
    integer(int32), parameter, public :: MODE_GROUP_ALL = CM_MODE_GROUP_ALL
        !! Mode code for taking all genes as a group in [[tox_gene_centroids(module):group_centroid(subroutine)]]

contains

    !> AUTHOR_LUKA_FAENSEN
    !| Computes the element-wise mean for a given set of vectors.
    pure subroutine mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, &
                                centroid, ierr)
        implicit none
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the input matrix.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), intent(in) :: n_selected_genes
            !! The number of genes in the current family to be averaged.
        integer(int32), dimension(n_selected_genes), intent(in) :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
        real(real64), dimension(n_axes), intent(out) :: centroid
            !! The output vector representing the computed centroid.
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Initialize error code
        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_genes, ierr, min=0_int32, max=n_genes, arg_pos=5_int32)

        if (is_err(ierr)) return

        call mean_vector_helper(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid)
    end subroutine mean_vector

    !> AUTHOR_LUKA_FAENSEN
    !| Computes the element-wise mean for a given set of vectors.
    pure subroutine mean_vector_helper(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid)
        implicit none
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the input matrix.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), intent(in) :: n_selected_genes
            !! The number of genes in the current family to be averaged.
        integer(int32), dimension(n_selected_genes), intent(in) :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
        real(real64), dimension(n_axes), intent(out) :: centroid
            !! The output vector representing the computed centroid.

        ! Local variables
        integer(int32) :: i_selected_gene, i_axis, gene_idx
        real(real64) :: sum_val

        ! If no genes are selected, return a zero vector
        centroid = 0.0_real64
        if (n_selected_genes == 0) return

        ! Compute the mean vector
        do concurrent (i_axis = 1:n_axes) local(sum_val) shared(gene_indices, centroid, n_selected_genes)
            sum_val = 0.0_real64
            ! For each selected gene, accumulate its expression value.
            do concurrent (i_selected_gene = 1:n_selected_genes) local(gene_idx) shared(gene_indices, expression_vectors) reduce(+:sum_val)
                gene_idx = gene_indices(i_selected_gene)
                sum_val = sum_val + expression_vectors(i_axis, gene_idx)
            end do
            ! Compute the mean for the current dimension by dividing through the number of selected genes.
            centroid(i_axis) = sum_val/real(n_selected_genes, real64)
        end do
    end subroutine mean_vector_helper

    !> AUTHOR_LUKA_FAENSEN
    !| Iterates over families, filters gene indices, and computes centroids.
    pure subroutine group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_group_indices, ierr, ortholog_set)
        implicit none
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! M_GENE_TO_FAM_DOC(expression_vectors)
        integer(int32), intent(in) :: mode
            !! used mode for grouping
            !!
            !! |      Method      |          Value            |
            !! |------------------|---------------------------|
            !! | Group Orthologs  |  CM_MODE_GROUP_ORTHOLOGS  |
            !! |    Group all     |     CM_MODE_GROUP_ALL     |
            !!
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! An output array for storing indices.
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).

        ! Initialize error code
        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=5_int32)

        call validate_all_in_range_int(gene_to_family, n_genes, ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL, arg_pos=4_int32)

        select case (mode)
            case (MODE_GROUP_ORTHOLOGS)
                if (.not. present(ortholog_set)) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=10_int32)
            case (MODE_GROUP_ALL)
            case default
                call set_err(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        end select

        if (is_err(ierr)) return

        call group_centroid_helper(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_group_indices, ortholog_set)
    end subroutine group_centroid

    !> AUTHOR_LUKA_FAENSEN
    !| Iterates over families, filters gene indices, and computes centroids.
    pure subroutine group_centroid_helper(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_group_indices, ortholog_set)
        implicit none
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! M_GENE_TO_FAM_DOC(expression_vectors)
        integer(int32), intent(in) :: mode
            !! used mode for grouping
            !!
            !! |      Method      |          Value            |
            !! |------------------|---------------------------|
            !! | Group Orthologs  |  CM_MODE_GROUP_ORTHOLOGS  |
            !! |    Group all     |     CM_MODE_GROUP_ALL     |
            !!
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! An output array for storing indices.
        logical, dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).

        ! Local variables
        integer(int32) :: i_gene, i_family, n_selected

        do i_family = 1, n_families
            ! Reset selected indices for the current family
            tmp_group_indices = 0
            n_selected = 0

            do i_gene = 1, n_genes
                ! Check if the gene belongs to the current family and in orthologs set if required
                if (gene_to_family(i_gene) == i_family) then
                    if (mode == MODE_GROUP_ORTHOLOGS) then
                        if (.not. ortholog_set(i_gene)) cycle
                    end if
                    n_selected = n_selected + 1
                    tmp_group_indices(n_selected) = i_gene
                end if
            end do

            call mean_vector_helper(expression_vectors, n_axes, n_genes, tmp_group_indices, n_selected, centroid_matrix(:, i_family))
        end do
    end subroutine group_centroid_helper

end module tox_gene_centroids

! =============================================================================
! C Wrapper Subroutine
! =============================================================================
!> C interface wrapper for mean_vector.
pure subroutine mean_vector_c(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr) &
    bind(c, name='mean_vector_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_gene_centroids, only: mean_vector
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_axes
        !! Number of axes (tissues/dimensions).
    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes in the input matrix.
    real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
        !! The input matrix of all gene expression vectors (n_axes x n_genes).
    integer(c_int), intent(in), target :: n_selected_genes
        !! The number of genes in the current family to be averaged.
    integer(c_int), dimension(n_selected_genes), intent(in), target :: gene_indices
        !! An array containing the column indices of the selected genes in 'expression_vectors'.
    real(c_double), dimension(n_axes), intent(out), target :: centroid_col
        !! The output vector representing the computed centroid.
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_selected_genes)
    M_CHECK_NON_NULL(expression_vectors)
    M_CHECK_NON_NULL(gene_indices)
    M_CHECK_NON_NULL(centroid_col)

    call mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr)
end subroutine mean_vector_c

!> C interface wrapper for group_centroid.
pure subroutine group_centroid_c(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, mode, ortholog_set, tmp_group_indices, ierr) &
    bind(c, name='group_centroid_c')
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_gene_centroids, only: group_centroid, MODE_GROUP_ORTHOLOGS, MODE_GROUP_ALL
    use tox_errors, only: is_err, set_err, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string, c_int_as_logical
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_axes
        !! Number of axes (tissues/dimensions).
    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes.
    integer(c_int), intent(in), target :: n_families
        !! Total number of families.
    real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
        !! The allocated length of the 'tmp_group_indices' array.
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
        !! M_GENE_TO_FAM_DOC(expression_vectors)
    character(c_char), dimension(9), intent(in), target :: mode
        !! A character array indicating the mode of operation ('orthologs' or 'all').
    real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
        !! Output matrix for centroids.
    integer(c_int), dimension(n_genes), intent(out), target :: tmp_group_indices
        !! Output array for selected indices.
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), dimension(n_genes), intent(in), target :: ortholog_set
        !! Integer array from C indicating subset membership.

    ! Local variables
    logical, dimension(:), allocatable :: ortholog_set_fortran
    integer(int32) :: mode_int
    character(len=:), allocatable :: mode_string

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(tmp_group_indices)
    M_CHECK_NON_NULL(expression_vectors)
    M_CHECK_NON_NULL(gene_to_family)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(centroid_matrix)
    M_CHECK_NON_NULL(tmp_group_indices)

    ! Convert raw character array to Fortran string
    call c_char_1d_as_string(mode, mode_string, ierr)
    if (is_err(ierr)) return

    ! If "orthologs" mode is selected, convert ortholog_set to logical
    ! If "all" mode is selected, call group_centroid directly without ortholog_set
    select case (mode_string)
    case ("orthologs")
        mode_int = MODE_GROUP_ORTHOLOGS
        M_CHECK_NON_NULL(ortholog_set)
        M_ALLOCATE(ortholog_set_fortran(n_genes))
        call c_int_as_logical(ortholog_set, ortholog_set_fortran)
    case ("all")
        mode_int = MODE_GROUP_ALL
    case default
        call set_err(ierr, ERR_INVALID_INPUT)
        return
    end select

    call group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroid_matrix, mode_int, tmp_group_indices, ierr, ortholog_set_fortran)
end subroutine group_centroid_c

