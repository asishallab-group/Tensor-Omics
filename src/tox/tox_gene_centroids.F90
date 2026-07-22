#include <src/macros.h>

!> Module for computing expression centroids of gene families.
!|
!| This module contains the core scientific kernel. The C and R interface
!| wrappers are defined outside the module for compatibility.
module tox_gene_centroids
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, validate_dimension_size, validate_in_range_int, validate_all_in_range_int
    use f42_utils, only: add_vector
    M_IMPLICIT_NONE

#define CM_MODE_GROUP_ORTHOLOGS 0_int32
#define CM_MODE_GROUP_ALL 1_int32

    integer(int32), parameter, public :: MODE_GROUP_ORTHOLOGS = CM_MODE_GROUP_ORTHOLOGS
        !! Mode code for grouping by orthologs in [[tox_gene_centroids(module):group_centroid(subroutine)]]
    integer(int32), parameter, public :: MODE_GROUP_ALL = CM_MODE_GROUP_ALL
        !! Mode code for taking all genes as a group in [[tox_gene_centroids(module):group_centroid(subroutine)]]

contains

    !> M_EXPORT_C
    !| summary: Computes the element-wise mean for a given set of vectors.
    !| AUTHOR_LUKA_FAENSEN
    pure subroutine mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, &
                                centroid, ierr)
        M_IMPLICIT_NONE
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
        call validate_all_in_range_int(gene_indices, n_selected_genes, ierr, min=1_int32, max=n_genes, arg_pos=4_int32)

        if (is_err(ierr)) return

        call mean_vector_helper(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid)
    end subroutine mean_vector

    !> AUTHOR_LUKA_FAENSEN
    !| Computes the element-wise mean for a given set of vectors.
    pure subroutine mean_vector_helper(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid)
        M_IMPLICIT_NONE
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
        do concurrent (i_selected_gene = 1:n_selected_genes) local(gene_idx) shared(gene_indices, centroid, expression_vectors)
            gene_idx = gene_indices(i_selected_gene)
            call add_vector(centroid, expression_vectors(:, gene_idx))
        end do

        do concurrent (i_axis = 1:n_axes) shared(centroid, n_selected_genes)
            centroid(i_axis) = centroid(i_axis) / real(n_selected_genes, real64)
        end do
    end subroutine mean_vector_helper

    !> M_EXPORT_C
    !| summary: Iterates over families, filters gene indices, and computes centroids.
    !| AUTHOR_LUKA_FAENSEN
    pure subroutine group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_group_indices, ierr, ortholog_set)
        M_IMPLICIT_NONE
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
            !! | Mode | Value |
            !! |------|-------|
            !! | Group Orthologs | [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]] |
            !! | Group all | [[tox_gene_centroids(module):MODE_GROUP_ALL(variable)]] |
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! An output array for storing indices.
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
            !! DM_REQUIRED_IF_MODE(mode, tox_gene_centroids, MODE_GROUP_ORTHOLOGS)

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
        M_IMPLICIT_NONE
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

        !TODO optimize: this rescans all n_genes for every family -> O(n_families * n_genes); a single
        !               bucketing pass (count genes per family, then fill per-family index lists) would be
        !               O(n_genes + n_families). Also, each family's centroid is independent work but the loop
        !               is serial and reuses a single shared `tmp_group_indices` scratch buffer across
        !               iterations, which would need to become per-iteration/private before this could be
        !               parallelized over families.
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

