#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shift_vectors(module)]]
!| The shift vector field: where each gene sits relative to its family's centroid.
!|
!| One vector per gene, from the centroid of the family it belongs to to the gene itself. It is
!| the input the relative-axis plane tools project and measure angles in, and what the
!| paralog-pattern detection reads a gene's direction and magnitude off.
module tox_shift_vectors_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: compute_shift_vector_field_c

contains

    !> summary: C-wrapper for [[tox_shift_vectors(module):compute_shift_vector_field(subroutine)]]
    !| Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
    subroutine compute_shift_vector_field_c(&
            n_tissues,&
            n_genes,&
            n_families,&
            expression_vectors,&
            family_centroids,&
            gene_to_fam,&
            shift_vectors,&
            ierr&
        ) bind(C, name="compute_shift_vector_field_c")
        use tox_shift_vectors, only: compute_shift_vector_field

        integer(c_int), intent(in), target :: n_tissues
            !! Expression vector dimension
        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of families
        real(c_double), dimension(n_tissues, n_genes), intent(in), target :: expression_vectors
            !! Gene expression matrix
        real(c_double), dimension(n_tissues, n_families), intent(in), target :: family_centroids
            !! Family centroid matrix
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_tissues, 2, n_genes), intent(out), target :: shift_vectors
            !! Output, real matrix array. For each gene it holds two vectors: the centroid of the gene's family first (a zero vector if no family is assigned), then the shift vector
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_tissues * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(shift_vectors, n_tissues * 2 * n_genes)

        call compute_shift_vector_field(&
            n_tissues = n_tissues,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            gene_to_fam = gene_to_fam,&
            shift_vectors = shift_vectors,&
            ierr = ierr&
        )
    end subroutine compute_shift_vector_field_c

end module tox_shift_vectors_c
#endif
