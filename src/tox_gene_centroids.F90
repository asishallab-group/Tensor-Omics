#include "macros.h"

!>Module for computing expression centroids of gene families.
!
! This module contains the core scientific kernel. The C and R interface
! wrappers are defined outside the module for compatibility.

module tox_gene_centroids
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once, is_ok
  implicit none

  integer(int32), parameter, public :: MODE_GROUP_ORTHOLOGS = 0
  integer(int32), parameter, public :: MODE_GROUP_ALL = 1

contains

  !> category: C-interface
  !| Computes the element-wise mean for a given set of vectors.
  pure subroutine mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes,&!asd
    centroid, ierr)
    implicit none
    !| Number of axes (tissues/dimensions).
    integer(int32), intent(in) :: n_axes
    !| Total number of genes in the input matrix.
    integer(int32), intent(in) :: n_genes
    !| The input matrix of all gene expression vectors (n_axes x n_genes).
    real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
    !| The number of genes in the current family to be averaged.
    integer(int32), intent(in) :: n_selected_genes
    !| An array containing the column indices of the selected genes in 'expression_vectors'.
    integer(int32), intent(in) :: gene_indices(n_selected_genes)
    !| The output vector representing the computed centroid.
    real(real64), intent(out) :: centroid(n_axes)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr

    ! Local variables
    integer(int32) :: i, j, gene_idx
    real(real64) :: sum_val

    ! Initialize error code
    call set_ok(ierr)

    ! Check for n_genes < 0
    if (n_axes <= 0 .or. n_genes <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if

    ! Check for invalid n_selected_genes (not < 0 and not > n_genes)
    if (n_selected_genes > n_genes) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! If no genes are selected, return a zero vector
    centroid = 0.0_real64
    if (n_selected_genes == 0) return

    ! Compute the mean vector
    do j = 1, n_axes
      sum_val = 0.0_real64
      ! For each selected gene, accumulate its expression value.
      do i = 1, n_selected_genes
        gene_idx = gene_indices(i)
        sum_val = sum_val + expression_vectors(j, gene_idx)
      end do
      ! Compute the mean for the current dimension by dividing through the number of selected genes.
      centroid(j) = sum_val / real(n_selected_genes, real64)
    end do
  end subroutine mean_vector

  !> category: C-interface
  !| Iterates over families, filters gene indices, and computes centroids.
  pure subroutine group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, mode, tmp_selected_indices, ierr, ortholog_set)
    implicit none
    integer(int32), intent(in) :: n_axes
      !! Number of axes (tissues/dimensions).
    integer(int32), intent(in) :: n_genes
      !! Total number of genes in the 'expression_vectors' matrix.
    integer(int32), intent(in) :: n_families
      !! Total number of gene families to compute centroids for.
    real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
      !! The input matrix of all gene expression vectors (n_axes x n_genes).
    integer(int32), intent(in) :: gene_to_family(n_genes)
      !! An array mapping each gene (by index) to a family ID.
    integer(int32), intent(in) :: mode
      !! used mode for grouping
      !!
      !! |       Mode       |                             Value                               |
      !! |------------------|-----------------------------------------------------------------|
      !! | Group Orthologs  |  [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]]  |
      !! |    Group all     |     [[tox_gene_centroids(module):MODE_GROUP_ALL(variable)]]     |
      !!
    real(real64), intent(out) :: centroid_matrix(n_axes, n_families)
      !! The output matrix (n_axes x n_families) to store the computed centroids.
    integer(int32), intent(out) :: tmp_selected_indices(n_genes)
      !! An output array for storing indices.
    integer(int32), intent(out) :: ierr
      !! Error code: 0 - success, non-zero = error
    logical, intent(in), optional :: ortholog_set(n_genes)
      !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
      !! DM_REQUIRED_IF_MODE(mode, tox_gene_centroids, MODE_GROUP_ORTHOLOGS)

    ! Local variables
    integer(int32) :: i, j, n_selected

    ! Initialize error code
    call set_ok(ierr)
    
    ! Determine the mode of operation
    if (mode /= MODE_GROUP_ALL .and. mode /= MODE_GROUP_ORTHOLOGS) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! If "orthologs" mode is selected, ensure ortholog_set is provided
    if (mode == MODE_GROUP_ORTHOLOGS .and. .not. present(ortholog_set)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! Check for arguments <= 0
    if (n_axes <= 0 .or. n_genes <= 0 .or. n_families <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if

    do j = 1, n_families
      ! Reset selected indices for the current family
      tmp_selected_indices = 0
      n_selected = 0

      do i = 1, n_genes
        ! Validate family ID
        if (gene_to_family(i) < 1 .or. gene_to_family(i) > n_families) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          return
        end if

        ! Check if the gene belongs to the current family and in orthologs set if required
        if (gene_to_family(i) == j) then
          if (mode == MODE_GROUP_ORTHOLOGS) then
            if (.not. ortholog_set(i)) cycle
          end if
          n_selected = n_selected + 1
          tmp_selected_indices(n_selected) = i
        end if
      end do

      call mean_vector(expression_vectors, n_axes, n_genes, tmp_selected_indices, n_selected, centroid_matrix(:, j), ierr)
      if (.not. is_ok(ierr)) return
    end do
  end subroutine group_centroid

end module tox_gene_centroids
