#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_gene_centroids(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_gene_centroids_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: mean_vector_c
    public :: group_centroid_orthologs_c
    public :: group_centroid_orthologs_expert_c
    public :: group_centroid_all_c
    public :: group_centroid_all_expert_c

contains

    !> summary: C-wrapper for [[tox_gene_centroids(module):mean_vector(subroutine)]]
    subroutine mean_vector_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_indices,&
            n_selected_genes,&
            centroid,&
            ierr&
        ) bind(C, name="mean_vector_c")
        use tox_gene_centroids, only: mean_vector

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the input matrix.
        integer(c_int), intent(in), target :: n_selected_genes
            !! The number of genes in the current family to be averaged.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_genes`.
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_selected_genes), intent(in), target :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        real(c_double), dimension(n_axes), intent(out), target :: centroid
            !! The output vector representing the computed centroid.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_selected_genes)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_indices, n_selected_genes)
        M_CHECK_ARRAY_NON_NULL(centroid, n_axes)

        call mean_vector(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_indices = gene_indices,&
            n_selected_genes = n_selected_genes,&
            centroid = centroid,&
            ierr = ierr&
        )
    end subroutine mean_vector_c

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid_orthologs(subroutine)]]
    subroutine group_centroid_orthologs_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            ortholog_set,&
            ierr&
        ) bind(C, name="group_centroid_orthologs_c")
        use tox_gene_centroids, only: group_centroid_orthologs

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families to compute centroids for.
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        logical(c_bool), dimension(n_genes), intent(in), target :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_genes) :: ortholog_set_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_genes)
        M_CHECK_ARRAY_NON_NULL(centroid_matrix, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(ortholog_set, n_genes)

        ortholog_set_f = ortholog_set

        call group_centroid_orthologs(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            ortholog_set = ortholog_set_f,&
            ierr = ierr&
        )
    end subroutine group_centroid_orthologs_c

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid_orthologs_expert(subroutine)]]
    subroutine group_centroid_orthologs_expert_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            tmp_group_indices,&
            ortholog_set,&
            ierr&
        ) bind(C, name="group_centroid_orthologs_expert_c")
        use tox_gene_centroids, only: group_centroid_orthologs_expert

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families to compute centroids for.
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_group_indices
            !! Work array for storing the indices of one family's genes.
        logical(c_bool), dimension(n_genes), intent(in), target :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_genes) :: ortholog_set_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_genes)
        M_CHECK_ARRAY_NON_NULL(centroid_matrix, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_group_indices, n_genes)
        M_CHECK_ARRAY_NON_NULL(ortholog_set, n_genes)

        ortholog_set_f = ortholog_set

        call group_centroid_orthologs_expert(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            tmp_group_indices = tmp_group_indices,&
            ortholog_set = ortholog_set_f,&
            ierr = ierr&
        )
    end subroutine group_centroid_orthologs_expert_c

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid_all(subroutine)]]
    subroutine group_centroid_all_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            ierr&
        ) bind(C, name="group_centroid_all_c")
        use tox_gene_centroids, only: group_centroid_all

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families to compute centroids for.
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_genes)
        M_CHECK_ARRAY_NON_NULL(centroid_matrix, n_axes * n_families)

        call group_centroid_all(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            ierr = ierr&
        )
    end subroutine group_centroid_all_c

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid_all_expert(subroutine)]]
    subroutine group_centroid_all_expert_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            tmp_group_indices,&
            ierr&
        ) bind(C, name="group_centroid_all_expert_c")
        use tox_gene_centroids, only: group_centroid_all_expert

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families to compute centroids for.
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_group_indices
            !! Work array for storing the indices of one family's genes.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_genes)
        M_CHECK_ARRAY_NON_NULL(centroid_matrix, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(tmp_group_indices, n_genes)

        call group_centroid_all_expert(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            tmp_group_indices = tmp_group_indices,&
            ierr = ierr&
        )
    end subroutine group_centroid_all_expert_c

end module tox_gene_centroids_c
#endif
