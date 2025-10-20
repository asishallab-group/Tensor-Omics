!>Module for computing expression centroids of gene families.
!
! This module contains the core scientific kernel. The C and R interface
! wrappers are defined outside the module for compatibility.

module tox_gene_centroids
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once, is_ok
  implicit none

  integer(int32), parameter, public :: GROUP_ORTHOLOGS = 0
  integer(int32), parameter, public :: GROUP_ALL = 1

contains

  !> Computes the element-wise mean for a given set of vectors.
  pure subroutine mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid, ierr)
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

  !> Iterates over families, filters gene indices, and computes centroids.
  pure subroutine group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, mode, selected_indices, ierr, ortholog_set)
    implicit none
    !| Number of axes (tissues/dimensions).
    integer(int32), intent(in) :: n_axes
    !| Total number of genes in the 'expression_vectors' matrix.
    integer(int32), intent(in) :: n_genes
    !| Total number of gene families to compute centroids for.
    integer(int32), intent(in) :: n_families
    !| The input matrix of all gene expression vectors (n_axes x n_genes).
    real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
    !| An array mapping each gene (by index) to a family ID.
    integer(int32), intent(in) :: gene_to_family(n_genes)
    !| A integer constant indicating the mode of operation ("all" or "orthologs").
    integer(int32), intent(in) :: mode
    !| The output matrix (n_axes x n_families) to store the computed centroids.
    real(real64), intent(out) :: centroid_matrix(n_axes, n_families)
    !| An output array for storing indices.
    integer(int32), intent(out) :: selected_indices(n_genes)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    !| A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
    logical, intent(in), optional :: ortholog_set(n_genes)

    ! Local variables
    integer(int32) :: i, j, n_selected

    ! Initialize error code
    call set_ok(ierr)
    
    ! Determine the mode of operation
    if (mode /= GROUP_ALL .and. mode /= GROUP_ORTHOLOGS) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! If "orthologs" mode is selected, ensure ortholog_set is provided
    if (mode == GROUP_ORTHOLOGS .and. .not. present(ortholog_set)) then
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
      selected_indices = 0
      n_selected = 0

      do i = 1, n_genes
        ! Validate family ID
        if (gene_to_family(i) < 1 .or. gene_to_family(i) > n_families) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          return
        end if

        ! Check if the gene belongs to the current family and in orthologs set if required
        if (gene_to_family(i) == j) then
          if (mode == GROUP_ORTHOLOGS) then
            if (.not. ortholog_set(i)) cycle
          end if
          n_selected = n_selected + 1
          selected_indices(n_selected) = i
        end if
      end do

      call mean_vector(expression_vectors, n_axes, n_genes, selected_indices, n_selected, centroid_matrix(:, j), ierr)
      if (.not. is_ok(ierr)) return
    end do
  end subroutine group_centroid

end module tox_gene_centroids

! =============================================================================
! C Wrapper Subroutine
! =============================================================================
!> C interface wrapper for mean_vector.
pure subroutine mean_vector_c(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr) &
  bind(c, name='mean_vector_c')
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_gene_centroids, only: mean_vector
  implicit none
  !| Number of axes (tissues/dimensions).
  integer(c_int), intent(in), value :: n_axes
  !| Total number of genes in the input matrix.
  integer(c_int), intent(in), value :: n_genes
  !| The input matrix of all gene expression vectors (n_axes x n_genes).
  real(c_double), intent(in) :: expression_vectors(n_axes, n_genes)
  !| The number of genes in the current family to be averaged.
  integer(c_int), intent(in), value :: n_selected_genes
  !| An array containing the column indices of the selected genes in 'expression_vectors'.
  integer(c_int), intent(in) :: gene_indices(n_selected_genes)
  !| The output vector representing the computed centroid.
  real(c_double), intent(out) :: centroid_col(n_axes)
  !| Error code: 0 - success, non-zero = error
  integer(c_int), intent(out) :: ierr

  call mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr)
end subroutine mean_vector_c

!> C interface wrapper for group_centroid.
pure subroutine group_centroid_c(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, mode, ortholog_set, selected_indices, &
                                 selected_indices_len, ierr) &
  bind(c, name='group_centroid_c')
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
  use tox_gene_centroids, only: group_centroid, GROUP_ORTHOLOGS, GROUP_ALL
  use tox_errors, only: is_ok, set_err, ERR_INVALID_INPUT
  use tox_conversions, only: c_char_1d_as_string, c_int_as_logical
  implicit none
  !| Number of axes (tissues/dimensions).
  integer(c_int), intent(in), value :: n_axes
  !| Total number of genes.
  integer(c_int), intent(in), value :: n_genes
  !| Total number of families.
  integer(c_int), intent(in), value :: n_families
  !| The allocated length of the 'selected_indices' array.
  integer(c_int), intent(in), value :: selected_indices_len
  !| Input expression vectors (passed from C).
  real(c_double), intent(in) :: expression_vectors(n_axes, n_genes)
  !| Array mapping gene index to family ID.
  integer(c_int), intent(in) :: gene_to_family(n_genes)
  !| A character array indicating the mode of operation ('orthologs' or 'all').
  character(c_char), intent(in) :: mode(10)
  !| Output matrix for centroids.
  real(c_double), intent(out) :: centroid_matrix(n_axes, n_families)
  !| Output array for selected indices.
  integer(c_int), intent(out) :: selected_indices(selected_indices_len)
  !| Error code: 0 - success, non-zero = error
  integer(c_int), intent(out) :: ierr
  !| Integer array from C indicating subset membership.
  integer(c_int), intent(in) :: ortholog_set(n_genes)

  ! Local variables
  logical :: ortholog_set_fortran(n_genes)
  integer :: i, mode_int
  character(len=:), allocatable :: mode_string

  ! Convert raw character array to Fortran string
  call c_char_1d_as_string(mode, mode_string, ierr)
  if (.not. is_ok(ierr)) return

  ! If "orthologs" mode is selected, convert ortholog_set to logical
  ! If "all" mode is selected, call group_centroid directly without ortholog_set
  select case (mode_string)
    case ("orthologs")
      mode_int = GROUP_ORTHOLOGS
      call c_int_as_logical(ortholog_set, ortholog_set_fortran)
    case ("all")
      mode_int = GROUP_ALL
    case default
      call set_err(ierr, ERR_INVALID_INPUT)
      return
  end select

  call group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                      centroid_matrix, mode_int, selected_indices, ierr, ortholog_set_fortran)
end subroutine group_centroid_c

! =============================================================================
! R Wrapper Subroutine
! =============================================================================
!> R interface wrapper for mean_vector.
pure subroutine mean_vector_r(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_gene_centroids, only: mean_vector
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
  real(real64), intent(out) :: centroid_col(n_axes)
  !| Error code: 0 - success, non-zero = error
  integer(int32), intent(out) :: ierr

  call mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr)
end subroutine mean_vector_r

!> R interface wrapper for group_centroid.
pure subroutine group_centroid_r(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, mode_raw, ortholog_set, &
                                 selected_indices, selected_indices_len, ierr)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_char
  use tox_gene_centroids, only: group_centroid, GROUP_ORTHOLOGS, GROUP_ALL
  use tox_conversions, only: c_char_1d_as_string
  use tox_errors, only: is_ok, set_err, ERR_INVALID_INPUT
  implicit none
  !| Number of axes (tissues/dimensions).
  integer(int32), intent(in) :: n_axes
  !| Total number of genes in the 'expression_vectors' matrix.
  integer(int32), intent(in) :: n_genes
  !| Total number of gene families to compute centroids for.
  integer(int32), intent(in) :: n_families
  !| The allocated length of the 'selected_indices' array.
  integer(int32), intent(in) :: selected_indices_len
  !| The input matrix of all gene expression vectors (n_axes x n_genes).
  real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
  !| An array mapping each gene (by index) to a family ID.
  integer(int32), intent(in) :: gene_to_family(n_genes)
  !| A raw character array indicating the mode of operation ("orthologs" or "all").
  character(c_char), intent(in) :: mode_raw(10)
  !| A logical array indicating if a gene is part of a specific subset.
  logical, intent(in) :: ortholog_set(n_genes)
  !| The output matrix (n_axes x n_families) to store the computed centroids.
  real(real64), intent(out) :: centroid_matrix(n_axes, n_families)
  !| An output array for storing selected gene indices.
  integer(int32), intent(out) :: selected_indices(selected_indices_len)
  !| Error code: 0 - success, non-zero = error
  integer(int32), intent(out) :: ierr

  ! Local variables
  integer(int32) :: mode_int
  character(len=:), allocatable :: mode_string

  ! Convert raw character array to Fortran string
  call c_char_1d_as_string(mode_raw, mode_string, ierr)
  if (.not. is_ok(ierr)) return

  ! Convert string to integer mode
  select case (mode_string)
    case ("orthologs")
      mode_int = GROUP_ORTHOLOGS
    case ("all")
      mode_int = GROUP_ALL
    case default
      call set_err(ierr, ERR_INVALID_INPUT)
      return
  end select
  call group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                      centroid_matrix, mode_int, selected_indices, ierr, ortholog_set)
end subroutine group_centroid_r
