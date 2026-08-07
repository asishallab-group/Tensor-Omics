#include <src/macros.h>

!> Implementation for computing the shift vector field for all genes.
!|
!| Hand-written implementation only. The generator turns this into the validating wrapper
!| [[tox_shift_vectors(module):compute_shift_vector_field]] in module `tox_shift_vectors`.
module tox_shift_vectors_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private
    public :: compute_shift_vector_field_impl
contains

    !> summary: Compute the shift vector field for all genes.
    !| AUTHOR_ALEXANDER_SCHWARZPAUL
    !| Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
    pure subroutine compute_shift_vector_field_impl(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                               gene_to_fam, shift_vectors)
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of families
        real(real64), intent(in) :: expression_vectors(n_tissues, n_genes)
            !! Gene expression matrix
        real(real64), intent(in) :: family_centroids(n_tissues, n_families)
            !! Family centroid matrix
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: shift_vectors(n_tissues, 2, n_genes)
            !! Output, real matrix array. For each gene it holds two vectors: the centroid of the gene's family first (a zero vector if no family is assigned), then the shift vector

        integer(int32) :: i_gene, fam_idx, i_tissue

        ! Each gene's shift vector only depends on its own expression and its family's centroid,
        ! so genes can be processed independently in any order.
        do concurrent (i_gene = 1:n_genes) local(fam_idx) shared(gene_to_fam, shift_vectors, family_centroids, expression_vectors, n_tissues)
            ! Get current centroid index from `i_gene`
            fam_idx = gene_to_fam(i_gene)

            if (fam_idx == M_GENE_TO_FAM_SENTINEL) then
                shift_vectors(:, 1, i_gene) = 0.0_real64
            else
                shift_vectors(:, 1, i_gene) = family_centroids(:, fam_idx)
            end if

            ! shift vector = gene's own expression minus its family centroid: it points from the
            ! family's average expression towards this gene, i.e. the direction and magnitude by
            ! which this paralog's expression has diverged from its family.
            do concurrent (i_tissue = 1:n_tissues) shared(shift_vectors, i_gene, expression_vectors, family_centroids, fam_idx)
                shift_vectors(i_tissue, 2, i_gene) = expression_vectors(i_tissue, i_gene) - family_centroids(i_tissue, fam_idx)
            end do
        end do
    end subroutine compute_shift_vector_field_impl
end module tox_shift_vectors_impl
