#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_gene_centroids(module)]]
!| Module for computing expression centroids of gene families.
!|
!| This module contains the core scientific kernel. The C and R bindings
!| wrappers are defined outside the module for compatibility.
module tox_gene_centroids_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: mean_vector_c
    public :: group_centroid_c

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
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(c_int), dimension(n_selected_genes), intent(in), target :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
        real(c_double), dimension(n_axes), intent(out), target :: centroid
            !! The output vector representing the computed centroid.
        integer(c_int), intent(out), target :: ierr
            !! Error code

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

    !> summary: C-wrapper for [[tox_gene_centroids(module):group_centroid(subroutine)]]
    subroutine group_centroid_c(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            mode,&
            tmp_group_indices,&
            ierr,&
            ortholog_set&
        ) bind(C, name="group_centroid_c")
        use tox_gene_centroids, only: group_centroid
        use tox_gene_centroids, only: MODE_GROUP_ALL, MODE_GROUP_ORTHOLOGS

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
        real(c_double), dimension(n_axes, n_families), intent(out), target :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        character(len=1, kind=c_char), dimension(15), intent(in), target :: mode
            !! used mode for grouping
            !!
            !! | Mode            | Value                                                         |
            !! |-----------------|---------------------------------------------------------------|
            !! | Group Orthologs | [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]] |
            !! | Group all       | [[tox_gene_centroids(module):MODE_GROUP_ALL(variable)]]       |
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_group_indices
            !! An output array for storing indices.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical(c_bool), dimension(n_genes), intent(in), optional :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
            !! This optional argument needs to be passed if used mode (`mode`) is [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]].
        integer(int32) :: mode_mode_f
        logical, dimension(:), allocatable :: ortholog_set_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_genes)
        M_CHECK_ARRAY_NON_NULL(centroid_matrix, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(mode, 15)
        M_CHECK_ARRAY_NON_NULL(tmp_group_indices, n_genes)

        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("group_orthologs")
                    mode_mode_f = MODE_GROUP_ORTHOLOGS
                case ("group_all")
                    mode_mode_f = MODE_GROUP_ALL
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        if (present(ortholog_set)) ortholog_set_f = ortholog_set

        call group_centroid(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            mode = mode_mode_f,&
            tmp_group_indices = tmp_group_indices,&
            ierr = ierr,&
            ortholog_set = ortholog_set_f&
        )
    end subroutine group_centroid_c

end module tox_gene_centroids_c
#endif
