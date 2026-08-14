#include <src/macros.h>

!> Expression centroids of gene families.
!|
!| `mean_vector` is the centroid of a set of expression vectors. `group_centroid_orthologs`
!| and `group_centroid_all` take the centroid of a family: over its orthologs only, or over
!| every gene in it -- two routines rather than one taking a flag, so which set a result is
!| over is visible at the call site.
!|
!| Generated from [[tox_gene_centroids_impl(module)]]; do not edit -- regenerate instead.
module tox_gene_centroids
    use tox_gene_centroids_impl, only: MODE_GROUP_ALL, MODE_GROUP_ORTHOLOGS, group_centroid_impl, mean_vector_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: mean_vector
    public :: group_centroid_orthologs
    public :: group_centroid_orthologs_expert
    public :: group_centroid_all
    public :: group_centroid_all_expert

contains

    !> summary: Validates its inputs, then calls [[tox_gene_centroids_impl(module):mean_vector_impl]].
    pure subroutine mean_vector(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_indices,&
            n_selected_genes,&
            centroid,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the input matrix.
        integer(int32), intent(in) :: n_selected_genes
            !! The number of genes in the current family to be averaged.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_genes`.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_selected_genes), intent(in) :: gene_indices
            !! An array containing the column indices of the selected genes in 'expression_vectors'.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        real(real64), dimension(n_axes), intent(out) :: centroid
            !! The output vector representing the computed centroid.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_genes, ierr, arg_pos=5_int32, min=0_int32, max=n_genes)
        call validate_all_in_range_real(expression_vectors, n_axes * n_genes, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(gene_indices, n_selected_genes, ierr, arg_pos=4_int32, min=1_int32, max=n_genes)
        if (is_err(ierr)) return
#endif

        call mean_vector_impl(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_indices = gene_indices,&
            n_selected_genes = n_selected_genes,&
            centroid = centroid&
        )
    end subroutine mean_vector

    !> summary: Validates its inputs, prepares what [[tox_gene_centroids_impl(module):group_centroid_impl]] needs, then calls it. The entry point to reach for first; see [[tox_gene_centroids(module):group_centroid_orthologs_expert]] to prepare it yourself.
    pure subroutine group_centroid_orthologs(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            ortholog_set,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        logical(c_bool), dimension(n_genes), intent(in) :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_group_indices

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_genes, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, arg_pos=4_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_group_indices(n_genes))

        call group_centroid_impl(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            mode = MODE_GROUP_ORTHOLOGS,&
            tmp_group_indices = tmp_group_indices,&
            ortholog_set = ortholog_set&
        )
    end subroutine group_centroid_orthologs

    !> summary: Validates its inputs, then calls [[tox_gene_centroids_impl(module):group_centroid_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_gene_centroids(module):group_centroid_orthologs]] does both.
    pure subroutine group_centroid_orthologs_expert(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            tmp_group_indices,&
            ortholog_set,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! Work array for storing the indices of one family's genes.
        logical(c_bool), dimension(n_genes), intent(in) :: ortholog_set
            !! A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_genes, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, arg_pos=4_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        call group_centroid_impl(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            mode = MODE_GROUP_ORTHOLOGS,&
            tmp_group_indices = tmp_group_indices,&
            ortholog_set = ortholog_set&
        )
    end subroutine group_centroid_orthologs_expert

    !> summary: Validates its inputs, prepares what [[tox_gene_centroids_impl(module):group_centroid_impl]] needs, then calls it. The entry point to reach for first; see [[tox_gene_centroids(module):group_centroid_all_expert]] to prepare it yourself.
    pure subroutine group_centroid_all(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_group_indices

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_genes, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, arg_pos=4_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_group_indices(n_genes))

        call group_centroid_impl(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            mode = MODE_GROUP_ALL,&
            tmp_group_indices = tmp_group_indices&
        )
    end subroutine group_centroid_all

    !> summary: Validates its inputs, then calls [[tox_gene_centroids_impl(module):group_centroid_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_gene_centroids(module):group_centroid_all]] does both.
    pure subroutine group_centroid_all_expert(&
            expression_vectors,&
            n_axes,&
            n_genes,&
            gene_to_family,&
            n_families,&
            centroid_matrix,&
            tmp_group_indices,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions).
        integer(int32), intent(in) :: n_genes
            !! Total number of genes in the 'expression_vectors' matrix.
        integer(int32), intent(in) :: n_families
            !! Total number of gene families to compute centroids for.
        real(real64), dimension(n_axes, n_genes), intent(in) :: expression_vectors
            !! The input matrix of all gene expression vectors (n_axes x n_genes).
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_axes, n_families), intent(out) :: centroid_matrix
            !! The output matrix (n_axes x n_families) to store the computed centroids.
        integer(int32), dimension(n_genes), intent(out) :: tmp_group_indices
            !! Work array for storing the indices of one family's genes.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_genes, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, arg_pos=4_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        call group_centroid_impl(&
            expression_vectors = expression_vectors,&
            n_axes = n_axes,&
            n_genes = n_genes,&
            gene_to_family = gene_to_family,&
            n_families = n_families,&
            centroid_matrix = centroid_matrix,&
            mode = MODE_GROUP_ALL,&
            tmp_group_indices = tmp_group_indices&
        )
    end subroutine group_centroid_all_expert

end module tox_gene_centroids
