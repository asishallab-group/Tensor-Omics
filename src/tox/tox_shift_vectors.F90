#include <src/macros.h>

!> Module for computing the shift vector field for all genes.
module tox_shift_vectors
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, is_err, validate_dimension_size, validate_all_in_range_int, validate_all_in_range_real
    M_IMPLICIT_NONE
contains

    !> M_EXPORT_C
    !| summary: Compute the shift vector field for all genes.
    !| AUTHOR_ALEXANDER_SCHWARZPAUL
    !| Computes the shift vectors by substracting the corresponding family centroid from the expression vector.
    pure subroutine compute_shift_vector_field(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                               gene_to_fam, shift_vectors, ierr)
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
        real(real64), intent(out) :: shift_vectors(n_tissues, 2, n_genes)
            !! Output, real matrix array, stores the centroid of the gene's family in `shift_vectors(:, 1, i_gene)` (zero vector if no family assigned) and the shift vectors in `shift_vectors(:, 2, i_gene)`
        integer(int32), intent(out) :: ierr
            !! Error code: 0 - success, non-zero = error

        integer(int32) :: i_gene, current_centroid, i

        call set_ok(ierr)

        call validate_dimension_size(n_tissues, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL, arg_pos=6_int32)
        call validate_all_in_range_real(expression_vectors, size(expression_vectors, kind=int32), ierr, arg_pos=4_int32)
        call validate_all_in_range_real(family_centroids, size(family_centroids, kind=int32), ierr, arg_pos=5_int32)

        if (is_err(ierr)) return

        call compute_shift_vector_field_helper(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                               gene_to_fam, shift_vectors)
    end subroutine compute_shift_vector_field

    !> AUTHOR_ALEXANDER_SCHWARZPAUL
    !| (no input validation) Compute the shift vector field for all genes.
    !| Computes the shift vectors by substracting the corresponding family centroid from the expression vector.
    pure subroutine compute_shift_vector_field_helper(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
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
        real(real64), intent(out) :: shift_vectors(n_tissues, 2, n_genes)
            !! Output, real matrix array, stores the centroid of the gene's family in `shift_vectors(:, 1, i_gene)` (zero vector if no family assigned) and the shift vectors in `shift_vectors(:, 2, i_gene)`

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
    end subroutine compute_shift_vector_field_helper
end module
