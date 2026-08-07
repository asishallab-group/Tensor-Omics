#include <src/macros.h>

!> Kernels for Euclidean distance computation for tensor omics.
!| The generator turns `euclidean_distance_impl` and `distance_to_centroid_impl` into the
!| validating wrappers `euclidean_distance` and `distance_to_centroid` in module
!| `tox_euclidean_distance`.
module tox_euclidean_distance_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE

#define DISTANCE_SENTINEL -1.0_real64

contains

    !> summary: Compute the Euclidean distance between two vectors.
    !| AUTHOR_VIVIAN_BASS
    !| Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
    pure subroutine euclidean_distance_impl(vec1, vec2, n_elements, result)
        integer(int32), intent(in) :: n_elements
            !! Dimension of both vectors
        real(real64), dimension(n_elements), intent(in) :: vec1
            !! First expression vector
        real(real64), dimension(n_elements), intent(in) :: vec2
            !! Second expression vector
        real(real64), intent(out) :: result
            !! Output scalar distance

        integer(int32) :: i_element
        real(real64) :: sum_squared_diff

        sum_squared_diff = 0.0_real64
        ! GFORTRAN BUG: do concurrent (i_element = 1:n_elements) shared(vec1, vec2) reduce(+:sum_squared_diff)
        do i_element = 1, n_elements
            sum_squared_diff = sum_squared_diff + (vec1(i_element) - vec2(i_element))**2
        end do
        result = sqrt(sum_squared_diff)
    end subroutine euclidean_distance_impl

    !> summary: Compute distance from each gene to its corresponding family centroid.
    !| AUTHOR_VIVIAN_BASS
    !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
    pure subroutine distance_to_centroid_impl(n_genes, n_families, genes, centroids, &
                                         gene_to_fam, distances, n_tissues)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        real(real64), dimension(n_tissues, n_genes), intent(in) :: genes
            !! Gene expression matrix (n_tissues × n_genes), column-major
        real(real64), dimension(n_tissues, n_families), intent(in) :: centroids
            !! Family centroid matrix (n_tissues × n_families), column-major, `DISTANCE_SENTINEL` for unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! M_GENE_TO_FAM_DOC(genes)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), dimension(n_genes), intent(out) :: distances
            !! Output distances array

        integer(int32) :: i_gene, family_idx

        do concurrent (i_gene = 1:n_genes) local(family_idx) shared(gene_to_fam, n_families, distances, genes, centroids, n_tissues)
            family_idx = gene_to_fam(i_gene)
            if (family_idx < 1 .or. family_idx > n_families) then
                distances(i_gene) = DISTANCE_SENTINEL  ! Error indicator
                cycle
            end if
            call euclidean_distance_impl(genes(:, i_gene), centroids(:, family_idx), n_tissues, distances(i_gene))
        end do
    end subroutine distance_to_centroid_impl
end module tox_euclidean_distance_impl
