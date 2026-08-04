#include <src/macros.h>

!> summary: Wrappers for [[tox_euclidean_distance_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_euclidean_distance
    use tox_euclidean_distance_kernel, only: distance_to_centroid_kernel, euclidean_distance_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, validate_all_in_range_int, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: euclidean_distance
    public :: distance_to_centroid

contains

    !> summary: Validates its inputs, then calls [[tox_euclidean_distance_kernel(module):euclidean_distance_kernel]].
    !| Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
    subroutine euclidean_distance(&
            vec1,&
            vec2,&
            n_elements,&
            result,&
            ierr&
        )
        integer(int32), intent(in) :: n_elements
            !! Dimension of both vectors
        real(real64), dimension(n_elements), intent(in) :: vec1
            !! First expression vector
        real(real64), dimension(n_elements), intent(in) :: vec2
            !! Second expression vector
        real(real64), intent(out) :: result
            !! Output scalar distance
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_elements, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(vec1, n_elements, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(vec2, n_elements, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call euclidean_distance_kernel(&
            vec1 = vec1,&
            vec2 = vec2,&
            n_elements = n_elements,&
            result = result&
        )
    end subroutine euclidean_distance

    !> summary: Validates its inputs, then calls [[tox_euclidean_distance_kernel(module):distance_to_centroid_kernel]].
    !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
    subroutine distance_to_centroid(&
            n_genes,&
            n_families,&
            genes,&
            centroids,&
            gene_to_fam,&
            distances,&
            n_tissues,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        real(real64), dimension(n_tissues, n_genes), intent(in) :: genes
            !! Gene expression matrix (n_tissues × n_genes), column-major
        real(real64), dimension(n_tissues, n_families), intent(in) :: centroids
            !! Family centroid matrix (n_tissues × n_families), column-major, `-1.0_real64` for unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_genes), intent(out) :: distances
            !! Output distances array
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(genes, n_tissues * n_genes, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(centroids, n_tissues * n_families, ierr, arg_pos=4_int32)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return

        call distance_to_centroid_kernel(&
            n_genes = n_genes,&
            n_families = n_families,&
            genes = genes,&
            centroids = centroids,&
            gene_to_fam = gene_to_fam,&
            distances = distances,&
            n_tissues = n_tissues&
        )
    end subroutine distance_to_centroid

end module tox_euclidean_distance
