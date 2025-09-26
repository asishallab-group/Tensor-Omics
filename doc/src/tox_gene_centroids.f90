!>Module for computing expression centroids of gene families.
!
! This module contains the core scientific kernel. The C and R interface
! wrappers are defined outside the module for compatibility.
module tox_gene_centroids
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once, is_ok
  implicit none
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
                                 centroid_matrix, use_all_mode, ortholog_set, selected_indices, ierr)
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
    !| The output matrix (n_axes x n_families) to store the computed centroids.
    real(real64), intent(out) :: centroid_matrix(n_axes, n_families)
    !| A logical flag; if true, all genes in a family are used.
    logical, intent(in) :: use_all_mode
    !| A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
    logical, intent(in) :: ortholog_set(n_genes)
    !| An output array for storing indices.
    integer(int32), intent(out) :: selected_indices(n_genes)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr

    ! Local variables
    integer(int32) :: i, j, n_selected
    integer(int32) :: local_selected_indices(n_genes)

    ! Initialize error code
    call set_ok(ierr)

    ! Check for n_genes < 0
    if (n_axes <= 0 .or. n_genes <= 0 .or. n_families <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if

    selected_indices = 0

    do j = 1, n_families
      n_selected = 0
      do i = 1, n_genes
        if (gene_to_family(i) < 1 .or. gene_to_family(i) > n_families) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          return
        end if
        if (gene_to_family(i) == j .and. (use_all_mode .or. ortholog_set(i))) then
          n_selected = n_selected + 1
          local_selected_indices(n_selected) = i
        end if
      end do
      call mean_vector(expression_vectors, n_axes, n_genes, local_selected_indices, n_selected, centroid_matrix(:, j), ierr)
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
  real(c_double), intent(in), target :: expression_vectors(n_axes, n_genes)
  !| The number of genes in the current family to be averaged.
  integer(c_int), intent(in), value :: n_selected_genes
  !| An array containing the column indices of the selected genes in 'expression_vectors'.
  integer(c_int), intent(in), target :: gene_indices(n_selected_genes)
  !| The output vector representing the computed centroid.
  real(c_double), intent(out), target :: centroid_col(n_axes)
  !| Error code: 0 - success, non-zero = error
  integer(c_int), intent(out) :: ierr

  call mean_vector(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid_col, ierr)
end subroutine mean_vector_c

!> C interface wrapper for group_centroid.
pure subroutine group_centroid_c(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                 centroid_matrix, use_all_mode, ortholog_set, &
                                 selected_indices, selected_indices_len, ierr) &
  bind(c, name='group_centroid_c')
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_gene_centroids, only: group_centroid
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
  real(c_double), intent(in), target :: expression_vectors(n_axes, n_genes)
  !| Array mapping gene index to family ID.
  integer(c_int), intent(in), target :: gene_to_family(n_genes)
  !| Integer flag from C (0=false, non-zero=true) to use all genes.
  integer(c_int), intent(in), value :: use_all_mode
  !| Integer array from C indicating subset membership.
  integer(c_int), intent(in), target :: ortholog_set(n_genes)
  !| Output matrix for centroids.
  real(c_double), intent(out), target :: centroid_matrix(n_axes, n_families)
  !| Output array for selected indices.
  integer(c_int), intent(out), target :: selected_indices(selected_indices_len)
  !| Error code: 0 - success, non-zero = error
  integer(c_int), intent(out) :: ierr

  ! Local variables
  logical :: use_all_mode_fortran
  logical :: ortholog_set_fortran(n_genes)
  integer :: i

  use_all_mode_fortran = (use_all_mode /= 0)
  do i = 1, n_genes
    ortholog_set_fortran(i) = (ortholog_set(i) /= 0)
  end do

  call group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                      centroid_matrix, use_all_mode_fortran, ortholog_set_fortran, selected_indices, ierr)
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
                                 centroid_matrix, use_all_mode, ortholog_set, &
                                 selected_indices, selected_indices_len, ierr)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_gene_centroids, only: group_centroid
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
  !| A logical flag; if true, all genes in a family are used.
  logical, intent(in) :: use_all_mode
  !| A logical array indicating if a gene is part of a specific subset.
  logical, intent(in) :: ortholog_set(n_genes)
  !| The output matrix (n_axes x n_families) to store the computed centroids.
  real(real64), intent(out) :: centroid_matrix(n_axes, n_families)
  !| An output array for storing selected gene indices.
  integer(int32), intent(out) :: selected_indices(selected_indices_len)
  !| Error code: 0 - success, non-zero = error
  integer(int32), intent(out) :: ierr

  call group_centroid(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                      centroid_matrix, use_all_mode, ortholog_set, selected_indices, ierr)
end subroutine group_centroid_r
