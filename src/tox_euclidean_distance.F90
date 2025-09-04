!> Module with Euclidean distance computation routines for tensor omics.
module tox_euclidean_distance
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

contains

  !> Compute the Euclidean distance between two vectors.
  !| Calculates the L2 norm: result = sqrt(sum((vec1_i - vec2_i)^2))
  pure subroutine euclidean_distance(vec1, vec2, d, result)
    !| Dimension of both vectors
    integer(int32), intent(in) :: d
    !| First expression vector
    real(real64), intent(in) :: vec1(d)
    !| Second expression vector
    real(real64), intent(in) :: vec2(d)
    !| Output scalar distance
    real(real64), intent(out) :: result
    integer(int32) :: i
    real(real64) :: sum_squared_diff
    sum_squared_diff = 0.0_real64
    do i = 1, d
      sum_squared_diff = sum_squared_diff + (vec1(i) - vec2(i))**2
    end do
    result = sqrt(sum_squared_diff)
  end subroutine euclidean_distance

  !> Compute distance from each gene to its corresponding family centroid.
  !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
  pure subroutine distance_to_centroid(n_genes, n_families, genes, centroids, &
                                       gene_to_fam, distances, d)
    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Total number of gene families
    integer(int32), intent(in) :: n_families
    !| Expression vector dimension
    integer(int32), intent(in) :: d
    !| Gene expression matrix (d × n_genes), column-major
    real(real64), intent(in) :: genes(d, n_genes)
    !| Family centroid matrix (d × n_families), column-major
    real(real64), intent(in) :: centroids(d, n_families)
    !| Gene-to-family mapping (1-based indexing)
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    !| Output distances array
    real(real64), intent(out) :: distances(n_genes)
    integer(int32) :: i, family_idx
    do i = 1, n_genes
      family_idx = gene_to_fam(i)
      if (family_idx < 1 .or. family_idx > n_families) then
        distances(i) = -1.0_real64  ! Error indicator
        cycle
      end if
      call euclidean_distance(genes(:, i), centroids(:, family_idx), d, distances(i))
    end do
  end subroutine distance_to_centroid




end module tox_euclidean_distance


!> R wrapper for euclidean_distance.
!| Calls euclidean_distance with standard Fortran types for R interface.
subroutine euclidean_distance_r(vec1, vec2, d, result)
  use tox_euclidean_distance
  !| Dimension of both vectors
  integer(int32), intent(in) :: d
  !| First expression vector
  real(real64), intent(in) :: vec1(d)
  !| Second expression vector
  real(real64), intent(in) :: vec2(d)
  !| Output scalar distance
  real(real64), intent(out) :: result
  call euclidean_distance(vec1, vec2, d, result)
end subroutine euclidean_distance_r

!> C wrapper for euclidean_distance.
!| Exposes euclidean_distance to C via iso_c_binding types.
subroutine euclidean_distance_c(vec1, vec2, d, result) bind(C, name="euclidean_distance_c")
  use iso_c_binding, only : c_int, c_double
  use tox_euclidean_distance
  !| Dimension of both vectors
  integer(c_int), intent(in), value :: d
  !| First expression vector
  real(c_double), intent(in), target :: vec1(d)
  !| Second expression vector
  real(c_double), intent(in), target :: vec2(d)
  !| Output scalar distance
  real(c_double), intent(out) :: result
  call euclidean_distance(vec1, vec2, d, result)
end subroutine euclidean_distance_c

!> R wrapper for distance_to_centroid.
!| Calls distance_to_centroid with standard Fortran types for R interface.
subroutine distance_to_centroid_r(n_genes, n_families, genes, centroids, &
                                  gene_to_fam, distances, d)
  use tox_euclidean_distance
  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Total number of gene families
  integer(int32), intent(in) :: n_families
  !| Expression vector dimension
  integer(int32), intent(in) :: d
  !| Gene expression matrix (d × n_genes), column-major
  real(real64), intent(in) :: genes(d, n_genes)
  !| Family centroid matrix (d × n_families), column-major
  real(real64), intent(in) :: centroids(d, n_families)
  !| Gene-to-family mapping (1-based indexing)
  integer(int32), intent(in) :: gene_to_fam(n_genes)
  !| Output distances array
  real(real64), intent(out) :: distances(n_genes)
  call distance_to_centroid(n_genes, n_families, genes, centroids, &
                                  gene_to_fam, distances, d)
end subroutine distance_to_centroid_r

!> C wrapper for distance_to_centroid.
!| Exposes distance_to_centroid to C via iso_c_binding types.
subroutine distance_to_centroid_c(n_genes, n_families, genes, centroids, & 
                                  gene_to_fam, distances, d) bind(C, name="distance_to_centroid_c")
  use iso_c_binding, only : c_int, c_double
  use tox_euclidean_distance
  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes
  !| Total number of gene families
  integer(c_int), intent(in), value :: n_families
  !| Expression vector dimension
  integer(c_int), intent(in), value :: d
  !| Gene expression matrix (d × n_genes), column-major
  real(c_double), intent(in), target :: genes(d, n_genes)
  !| Family centroid matrix (d × n_families), column-major
  real(c_double), intent(in), target :: centroids(d, n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  !| Output distances array
  real(c_double), intent(out), target :: distances(n_genes)
  
  call distance_to_centroid(n_genes, n_families, genes, centroids, &
                            gene_to_fam, distances, d)
end subroutine distance_to_centroid_c