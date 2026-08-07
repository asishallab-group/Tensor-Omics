#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_euclidean_distance(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_euclidean_distance_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: euclidean_distance_c
    public :: distance_to_centroid_c

contains

    !> summary: C-wrapper for [[tox_euclidean_distance(module):euclidean_distance(subroutine)]]
    !| Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
    subroutine euclidean_distance_c(&
            vec1,&
            vec2,&
            n_elements,&
            result,&
            ierr&
        ) bind(C, name="euclidean_distance_c")
        use tox_euclidean_distance, only: euclidean_distance

        integer(c_int), intent(in), target :: n_elements
            !! Dimension of both vectors
        real(c_double), dimension(n_elements), intent(in), target :: vec1
            !! First expression vector
        real(c_double), dimension(n_elements), intent(in), target :: vec2
            !! Second expression vector
        real(c_double), intent(out), target :: result
            !! Output scalar distance
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_elements)
        M_CHECK_NON_NULL(result)
        M_CHECK_ARRAY_NON_NULL(vec1, n_elements)
        M_CHECK_ARRAY_NON_NULL(vec2, n_elements)

        call euclidean_distance(&
            vec1 = vec1,&
            vec2 = vec2,&
            n_elements = n_elements,&
            result = result,&
            ierr = ierr&
        )
    end subroutine euclidean_distance_c

    !> summary: C-wrapper for [[tox_euclidean_distance(module):distance_to_centroid(subroutine)]]
    !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
    subroutine distance_to_centroid_c(&
            n_genes,&
            n_families,&
            genes,&
            centroids,&
            gene_to_fam,&
            distances,&
            n_tissues,&
            ierr&
        ) bind(C, name="distance_to_centroid_c")
        use tox_euclidean_distance, only: distance_to_centroid

        integer(c_int), intent(in), target :: n_genes
            !! Total number of genes
        integer(c_int), intent(in), target :: n_families
            !! Total number of gene families
        integer(c_int), intent(in), target :: n_tissues
            !! Expression vector dimension
        real(c_double), dimension(n_tissues, n_genes), intent(in), target :: genes
            !! Gene expression matrix (n_tissues × n_genes), column-major
        real(c_double), dimension(n_tissues, n_families), intent(in), target :: centroids
            !! Family centroid matrix (n_tissues × n_families), column-major, `-1.0_real64` for unassigned genes
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_genes), intent(out), target :: distances
            !! Output distances array
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_ARRAY_NON_NULL(genes, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(centroids, n_tissues * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(distances, n_genes)

        call distance_to_centroid(&
            n_genes = n_genes,&
            n_families = n_families,&
            genes = genes,&
            centroids = centroids,&
            gene_to_fam = gene_to_fam,&
            distances = distances,&
            n_tissues = n_tissues,&
            ierr = ierr&
        )
    end subroutine distance_to_centroid_c

end module tox_euclidean_distance_c
#endif
