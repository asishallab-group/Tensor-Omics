#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_gene_centroids(module)]]
!| Module for computing expression centroids of gene families.
module tox_gene_centroids_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[tox_gene_centroids(module):mean_vector(subroutine)]]
    !| Computes the element-wise mean for a given set of vectors.
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
        use tox_gene_centroids
        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the input matrix.
        integer(c_int), intent(in), target :: n_selected_genes
            !! The number of genes in the current family to be averaged.
        real(c_double), intent(in), dimension(n_axes, n_genes), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), intent(in), dimension(n_selected_genes), target :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
        real(c_double), intent(out), dimension(n_axes), target :: centroid
            !! The output vector representing the computed centroid.
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0 - success, non-zero = error
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(expression_vectors)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(gene_indices)
        M_CHECK_NON_NULL(n_selected_genes)
        M_CHECK_NON_NULL(centroid)
        call  mean_vector(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_indices = gene_indices,&
            n_selected_genes = n_selected_genes,&
            centroid = centroid,&
            ierr = ierr&
        )
    end subroutine mean_vector_c

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid(subroutine)]]
    !| Iterates over families, filters gene indices, and computes centroids.
    subroutine group_centroid_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            mode,&
            tmp_selected_indices,&
            ierr,&
            ortholog_set&
            ) bind(C, name="group_centroid_c")
        use tox_gene_centroids, only: group_centroid
        use tox_gene_centroids
        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families to compute centroids for.
        real(c_double), intent(in), dimension(n_axes, n_genes), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), intent(in), dimension(n_genes), target :: gene_to_family
            !! An array mapping each gene (by index) to a family ID.
        real(c_double), intent(out), dimension(n_axes, n_families), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        character(len=1, kind=c_char), intent(in), dimension(15), target :: mode
            !! used mode for grouping
            !! 
            !! |       Mode       |                             Value                               |
            !! |------------------|-----------------------------------------------------------------|
            !! | Group Orthologs  |   "group_orthologs"    |
            !! |    Group all     |      "group_all"       |
        integer(c_int), intent(out), dimension(n_genes), target :: tmp_selected_indices
            !! An output array for storing indices.
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0 - success, non-zero = error
        integer(c_int), intent(in), dimension(n_genes), target :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
            !! This optional argument needs to be passed if used mode (`mode`) is [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]].
        character(len=:), allocatable :: mode_f
        logical, allocatable, dimension(:) :: ortholog_set_f
        integer(int32) :: mode_int_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(expression_vectors)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(gene_to_family)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(centroid_matrix)
        M_CHECK_NON_NULL(mode)
        M_CHECK_NON_NULL(tmp_selected_indices)
        call c_char_1d_as_string(mode, mode_f, ierr)
        if (is_err(ierr)) return
        select case (mode_f)
            case ("group_orthologs")
                    mode_int_f = MODE_GROUP_ORTHOLOGS
            case ("group_all")
                    mode_int_f = MODE_GROUP_ALL
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
                return
        end select
        if (mode_int_f == MODE_GROUP_ORTHOLOGS) then
            M_CHECK_NON_NULL(ortholog_set)
            M_ALLOCATE(ortholog_set_f(n_genes))
            call c_int_as_logical(ortholog_set, ortholog_set_f)
        end if
        if (mode_int_f == MODE_GROUP_ORTHOLOGS) then
            call  group_centroid(&
                expression_vectors = expression_vectors,&
                n_axes = n_axes,&
                n_genes = n_genes,&
                gene_to_family = gene_to_family,&
                n_families = n_families,&
                centroid_matrix = centroid_matrix,&
                mode = mode_int_f,&
                tmp_selected_indices = tmp_selected_indices,&
                ierr = ierr,&
                ortholog_set = ortholog_set_f&
            )
        else
            call  group_centroid(&
                expression_vectors = expression_vectors,&
                n_axes = n_axes,&
                n_genes = n_genes,&
                gene_to_family = gene_to_family,&
                n_families = n_families,&
                centroid_matrix = centroid_matrix,&
                mode = mode_int_f,&
                tmp_selected_indices = tmp_selected_indices,&
                ierr = ierr&
            )
        end if
    end subroutine group_centroid_c

end module tox_gene_centroids_c
#endif