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
    use f42_vector_impl, only: add_vector
    M_IMPLICIT_NONE

    private
    public :: mean_vector_impl, group_centroid_impl
    public :: MODE_GROUP_ORTHOLOGS, MODE_GROUP_ALL

    integer(int32), parameter :: MODE_GROUP_ORTHOLOGS = 0_int32
        !! Mode code for grouping by orthologs
    integer(int32), parameter :: MODE_GROUP_ALL = 1_int32
        !! Mode code for taking all genes as a group

contains

    !> summary: Computes the element-wise mean for a given set of vectors.
    !| AUTHOR_LUKA_FAENSEN
    pure subroutine mean_vector_impl(expression_vectors, n_axes, n_genes, gene_indices, n_selected_genes, centroid)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the input matrix.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), intent(in) :: n_selected_genes
            !! The number of genes in the current family to be averaged.
            !! DM_MIN(0_int32)
            !! DM_MAX(n_genes)
        integer(int32), dimension(n_selected_genes), intent(in) :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes)
        real(real64), dimension(n_axes), intent(out) :: centroid
            !! The output vector representing the computed centroid.

        integer(int32) :: i_selected_gene, i_axis, gene_idx

        ! If no genes are selected, return a zero vector
        centroid = 0.0_real64
        if (n_selected_genes == 0) return

        do concurrent (i_selected_gene = 1:n_selected_genes) local(gene_idx) shared(gene_indices, centroid, expression_vectors)
            gene_idx = gene_indices(i_selected_gene)
            call add_vector(centroid, expression_vectors(:, gene_idx))
        end do

        do concurrent (i_axis = 1:n_axes) shared(centroid, n_selected_genes)
            centroid(i_axis) = centroid(i_axis) / real(n_selected_genes, real64)
        end do
    end subroutine mean_vector_impl

    !> summary: Iterates over families, filters gene indices, and computes centroids.
    !| AUTHOR_LUKA_FAENSEN
    pure subroutine group_centroid_impl(expression_vectors, n_axes, n_genes, gene_to_family, n_families, &
                                   centroid_matrix, mode, tmp_group_indices, ortholog_set)
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
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! Work array for storing the indices of one family's genes.
        logical(c_bool), dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
            !! DM_REQUIRED_IF_MODE(mode, tox_gene_centroids_impl, MODE_GROUP_ORTHOLOGS)

        integer(int32) :: i_gene, i_family, n_selected

        do i_family = 1, n_families
            ! Reset selected indices for the current family
            tmp_group_indices = 0
            n_selected = 0

            do i_gene = 1, n_genes
                ! Check if the gene belongs to the current family and, if required, is an ortholog
                if (gene_to_family(i_gene) == i_family) then
                    if (mode == MODE_GROUP_ORTHOLOGS) then
                        if (.not. ortholog_set(i_gene)) cycle
                    end if
                    n_selected = n_selected + 1
                    tmp_group_indices(n_selected) = i_gene
                end if
            end do

            call mean_vector_impl(expression_vectors, n_axes, n_genes, tmp_group_indices, n_selected, centroid_matrix(:, i_family))
        end do
    end subroutine group_centroid_impl

end module tox_gene_centroids_impl
