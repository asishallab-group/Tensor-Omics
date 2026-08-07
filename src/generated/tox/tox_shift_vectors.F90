#include <src/macros.h>

!> summary: Wrappers for [[tox_shift_vectors_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_shift_vectors
    use tox_shift_vectors_impl, only: compute_shift_vector_field_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, validate_all_in_range_int, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: compute_shift_vector_field

contains

    !> summary: Validates its inputs, then calls [[tox_shift_vectors_impl(module):compute_shift_vector_field_impl]].
    !| Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
    subroutine compute_shift_vector_field(&
            n_tissues,&
            n_genes,&
            n_families,&
            expression_vectors,&
            family_centroids,&
            gene_to_fam,&
            shift_vectors,&
            ierr&
        )
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of families
        real(real64), dimension(n_tissues, n_genes), intent(in) :: expression_vectors
            !! Gene expression matrix
        real(real64), dimension(n_tissues, n_families), intent(in) :: family_centroids
            !! Family centroid matrix
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_tissues, 2, n_genes), intent(out) :: shift_vectors
            !! Output, real matrix array. For each gene it holds two vectors: the centroid of the gene's family first (a zero vector if no family is assigned), then the shift vector
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_tissues, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(expression_vectors, n_tissues * n_genes, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(family_centroids, n_tissues * n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=6_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        call compute_shift_vector_field_impl(&
            n_tissues = n_tissues,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            gene_to_fam = gene_to_fam,&
            shift_vectors = shift_vectors&
        )
    end subroutine compute_shift_vector_field

end module tox_shift_vectors
