!> Module for computing the shift vector field for all genes.
module tox_shift_vectors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once
contains

  !> Compute the shift vector field for all genes.
  !| Computes the shift vectors by substracting the corresponding family centroid from the expression vector.
  pure subroutine compute_shift_vector_field(d, n_genes, n_families, expression_vectors, family_centroids, &
                                             gene_to_centroid, shift_vectors, ierr)
    implicit none

    !| Expression vector dimension
    integer(int32), intent(in) :: d
    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Total number of families
    integer(int32), intent(in) :: n_families
    !| Gene expression matrix (d × n_genes)
    real(real64), intent(in) :: expression_vectors(d, n_genes)
    !| Family centroid matrix (d × n_families)
    real(real64), intent(in) :: family_centroids(d, n_families)
    !| Mapping from genes to family centroids to genes of `expression_vectors` -> `family_centroids(:, gene_to_centroid(n))` returns the family centroid of `expression_vectors(:,n)`
    integer(int32), intent(in) :: gene_to_centroid(n_genes)
    !| Output, real matrix array, size = 2d x `n_genes`, stores the centroid of the gene's family in rows 1..d and the shift vectors in rows d+1...2d
    real(real64), intent(out) :: shift_vectors(2*d, n_genes)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr

    !| Local variables
    integer(int32) :: current_gene, current_centroid, i

    ! Initialize error code
    call set_ok(ierr)

    ! Check for correct 0 dimension
    if (d == 0 .or. n_genes == 0 .or. n_families == 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if

    ! For each gene do
    do current_gene = 1, n_genes
      ! Check if `gene_to_centroid` mapping is in valid range (length of `family_centroids` array)
      if (gene_to_centroid(current_gene) < 1 .or. gene_to_centroid(current_gene) > n_families) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
      end if
      ! Get current centroid index from `current_gene`
      current_centroid = gene_to_centroid(current_gene)
      ! Copy family centroid to first half of the `shift_vectors`
      shift_vectors(1:d, current_gene) = family_centroids(:, current_centroid)
      ! Substract family centroid from `expression_vector` and write it to second half of the `shift_vectors`
      do i = 1, d
        shift_vectors(d + i, current_gene) = expression_vectors(i, current_gene) - family_centroids(i, current_centroid)
      end do
    end do
  end subroutine
end module

!> R wrapper for `compute_shift_vector_field`
!| Calls `compute_shift_vector_field` with standard Fortran types for R interface.
!| When using these R wrapper functions, copies of the arrays will be created. No direct modification of the original R objects occurs..
pure subroutine compute_shift_vector_field_r(d, n_genes, n_families, expression_vectors, family_centroids, &
                                             gene_to_centroid, shift_vectors, ierr)
  use tox_shift_vectors

  !| Expression vector dimension
  integer(int32), intent(in) :: d
  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Total number of families
  integer(int32), intent(in) :: n_families
  !| Gene expression matrix `(d × n_genes)`
  real(real64), intent(in) :: expression_vectors(d, n_genes)
  !| Family centroid matrix `(d × n_families)`
  real(real64), intent(in) :: family_centroids(d, n_families)
  !| Mapping from genes to family centroids to genes of `expression_vectors` -> `family_centroids(:, gene_to_centroid(n))` returns the family centroid of `expression_vectors(:,n)`
  integer(int32), intent(in) :: gene_to_centroid(n_genes)
  !| Output, real matrix array, size = 2d x `n_genes`, stores the centroid of the gene's family in rows 1..d and the shift vectors in rows d+1...2d
  real(real64), intent(out) :: shift_vectors(2*d, n_genes)
  !| Error code: 0 - success, non-zero = error
  integer(int32), intent(out) :: ierr
  call compute_shift_vector_field(d, n_genes, n_families, expression_vectors, family_centroids, &
                                  gene_to_centroid, shift_vectors, ierr)
end subroutine compute_shift_vector_field_r

!> C wrapper for `compute_shift_vector_field`.
!| Exposes `compute_shift_vector_field` to C via iso_c_binding types with explicit dimensions.
!| When using these C wrapper functions, no copies of the arrays will be created. The Fortran routine will operate directly on the memory provided by the caller.
pure subroutine compute_shift_vector_field_c(d, n_genes, n_families, expression_vectors, family_centroids, gene_to_centroid, &
                                             shift_vectors, ierr) bind(C, name="compute_shift_vector_field_c")
  use iso_c_binding
  use tox_shift_vectors
  !| Expression vector dimension
  integer(c_int), intent(in), value :: d
  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes
  !| Total number of families
  integer(c_int), intent(in), value :: n_families
  !| Gene expression matrix `(d × n_genes)`
  real(c_double), intent(in), target :: expression_vectors(d, n_genes)
  !| Family centroid matrix `(d × n_families)`
  real(c_double), intent(in), target :: family_centroids(d, n_families)
  !| Mapping from genes to families (family IDs for each gene in `expression_vectors`)
  integer(c_int), intent(in), target :: gene_to_centroid(n_genes)
  !| Output, real matrix array, size = 2d x `n_genes`, stores the centroid of the gene's family in rows 1..d and the shift vectors in rows d+1...2d
  real(c_double), intent(out), target :: shift_vectors(2*d, n_genes)
  !| Error code: 0 - success, non-zero = error
  integer(c_int), intent(out) :: ierr

  call compute_shift_vector_field(d, n_genes, n_families, expression_vectors, family_centroids, &
                                  gene_to_centroid, shift_vectors, ierr)
end subroutine compute_shift_vector_field_c
