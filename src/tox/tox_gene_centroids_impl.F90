#include <src/macros.h>

!> Expression centroids of gene families.
!|
!| `mean_vector` is the centroid of a set of expression vectors. `group_centroid_orthologs`
!| and `group_centroid_all` take the centroid of a family: over its orthologs only, or over
!| every gene in it -- two routines rather than one taking a flag, so which set a result is
!| over is visible at the call site.
module tox_gene_centroids_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    M_IMPLICIT_NONE

    private
    public :: mean_vector_impl, group_centroid_impl
    public :: MODE_GROUP_ORTHOLOGS, MODE_GROUP_ALL

    integer(int32), parameter :: MODE_GROUP_ORTHOLOGS = 0_int32
        !! Mode code for grouping by orthologs
    integer(int32), parameter :: MODE_GROUP_ALL = 1_int32
        !! Mode code for taking all genes as a group

contains

    !> summary: Computes the element-wise mean of the selected gene vectors.
    !| AUTHOR_LUKA_FAENSEN
    !| A selection without genes gives the zero vector.
    pure subroutine mean_vector_impl(expression_vectors, n_axes, n_genes, genes_selection_mask, n_selected_genes, centroid)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the input matrix.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        logical(c_bool), dimension(n_genes), intent(in) :: genes_selection_mask
            !! `.true.` for the genes (columns of `expression_vectors`) to average
        integer(int32), intent(in) :: n_selected_genes
            !! count of `.true.` values in `genes_selection_mask`
            !! DM_MIN(0_int32)
        real(real64), dimension(n_axes), intent(out) :: centroid
            !! The output vector representing the computed centroid.

        integer(int32) :: i_gene, i_axis, scale_exponent
        real(real64) :: scaled_count

        ! If no genes are selected, return a zero vector
        centroid = 0.0_real64
        if (n_selected_genes == 0) return

        ! A mean is never larger than its largest value, but the plain sum can overflow where the
        ! mean does not. So every value is scaled by 2**-scale_exponent before it is added, where
        ! 2**scale_exponent > n_selected_genes: exact, and the scaled sum stays below the largest
        ! value. Dividing by the equally scaled count then gives the same quotient, with the same
        ! rounding, as the plain sum divided by the count.
        ! The genes are summed in order, not in a do concurrent: every gene adds into the same
        ! centroid, which concurrent iterations may not do. The axes are independent.
        scale_exponent = exponent(real(n_selected_genes, real64))
        do i_gene = 1, n_genes
            if (.not. genes_selection_mask(i_gene)) cycle
            do concurrent (i_axis = 1:n_axes) shared(centroid, expression_vectors, i_gene, scale_exponent)
                centroid(i_axis) = centroid(i_axis) + scale(expression_vectors(i_axis, i_gene), -scale_exponent)
            end do
        end do

        scaled_count = scale(real(n_selected_genes, real64), -scale_exponent)
        do concurrent (i_axis = 1:n_axes) shared(centroid, scaled_count)
            centroid(i_axis) = centroid(i_axis)/scaled_count
        end do
    end subroutine mean_vector_impl

    !> summary: Computes one centroid per gene family, of all its genes or of its orthologs only.
    !| AUTHOR_LUKA_FAENSEN
    !| A family without selected genes gets the zero vector.
    pure subroutine group_centroid_impl(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_family_genes, ortholog_set)
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
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        integer(int32), intent(in) :: mode
            !! used mode for grouping
            !!
            !! | Mode | Value | Procedure |
            !! |------|-------|-----------|
            !! | Group orthologs | [[tox_gene_centroids_impl(module):MODE_GROUP_ORTHOLOGS(variable)]] | group_centroid_orthologs |
            !! | Group all | [[tox_gene_centroids_impl(module):MODE_GROUP_ALL(variable)]] | group_centroid_all |
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        logical(c_bool), dimension(n_genes), intent(out) :: tmp_family_genes
            !! Work array: `.true.` for the genes of the family being averaged.
        logical(c_bool), dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
            !! DM_REQUIRED_IF_MODE(mode, tox_gene_centroids_impl, MODE_GROUP_ORTHOLOGS)

        integer(int32) :: i_gene, i_family, n_selected

        do i_family = 1, n_families
            ! select the genes of the current family and, if required, only its orthologs
            n_selected = 0
            do i_gene = 1, n_genes
                tmp_family_genes(i_gene) = gene_to_family(i_gene) == i_family
                if (mode == MODE_GROUP_ORTHOLOGS .and. tmp_family_genes(i_gene)) then
                    tmp_family_genes(i_gene) = ortholog_set(i_gene)
                end if
                if (tmp_family_genes(i_gene)) n_selected = n_selected + 1
            end do

            call mean_vector_impl(expression_vectors, n_axes, n_genes, tmp_family_genes, n_selected, centroid_matrix(:, i_family))
        end do
    end subroutine group_centroid_impl

end module tox_gene_centroids_impl
