#include <src/macros.h>

!> Module with Euclidean distance computation routines for tensor omics.
module tox_euclidean_distance
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: is_err, set_err, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, set_ok
    implicit none

#define DISTANCE_SENTINEL -1.0_real64

contains

    !> AUTHOR_VIVIAN_BASS
    !| Compute the Euclidean distance between two vectors.
    !| Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
    pure subroutine euclidean_distance(vec1, vec2, n_elements, result, ierr)
        integer(int32), intent(in) :: n_elements
            !! Dimension of both vectors
        real(real64), dimension(n_elements), intent(in) :: vec1
            !! First expression vector
        real(real64), dimension(n_elements), intent(in) :: vec2
            !! Second expression vector
        real(real64), intent(out) :: result
            !! Output scalar distance
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_elements, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(vec1, n_elements, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(vec2, n_elements, ierr, arg_pos=2_int32)

        if (is_err(ierr)) return

        call euclidean_distance_helper(vec1, vec2, n_elements, result)
    end subroutine euclidean_distance

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Compute the Euclidean distance between two vectors.
    !| Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
    pure subroutine euclidean_distance_helper(vec1, vec2, n_elements, result)
        integer(int32), intent(in) :: n_elements
            !! Dimension of both vectors
        real(real64), dimension(n_elements), intent(in) :: vec1
            !! First expression vector
        real(real64), dimension(n_elements), intent(in) :: vec2
            !! Second expression vector
        real(real64), intent(out) :: result
            !! Output scalar distance

        integer(int32) :: i_element
        real(real64) :: sum_squared_diff

        sum_squared_diff = 0.0_real64
        ! GFORTRAN BUG: do concurrent (i_element = 1:n_elements) shared(vec1, vec2) reduce(+:sum_squared_diff)
        do i_element = 1, n_elements
            sum_squared_diff = sum_squared_diff + (vec1(i_element) - vec2(i_element))**2
        end do
        result = sqrt(sum_squared_diff)
    end subroutine euclidean_distance_helper

    !> AUTHOR_VIVIAN_BASS
    !| Compute distance from each gene to its corresponding family centroid.
    !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
    pure subroutine distance_to_centroid(n_genes, n_families, genes, centroids, &
                                         gene_to_fam, distances, n_tissues, ierr)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        real(real64), dimension(n_tissues, n_genes), intent(in) :: genes
            !! Gene expression matrix (n_tissues × n_genes), column-major
        real(real64), dimension(n_tissues, n_families), intent(in) :: centroids
            !! Family centroid matrix (n_tissues × n_families), column-major, `DISTANCE_SENTINEL` for unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! M_GENE_TO_FAM_DOC(genes)
        real(real64), dimension(n_genes), intent(out) :: distances
            !! Output distances array
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_tissues, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(genes, size(genes, kind=int32), ierr, arg_pos=3_int32)
        call validate_all_in_range_real(centroids, size(centroids, kind=int32), ierr, arg_pos=4_int32)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL, arg_pos=5_int32)

        if (is_err(ierr)) return

        call distance_to_centroid_helper(n_genes, n_families, genes, centroids, gene_to_fam, distances, n_tissues)
    end subroutine distance_to_centroid

    !> AUTHOR_VIVIAN_BASS
    !| Compute distance from each gene to its corresponding family centroid.
    !| For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
    pure subroutine distance_to_centroid_helper(n_genes, n_families, genes, centroids, &
                                         gene_to_fam, distances, n_tissues)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        real(real64), dimension(n_tissues, n_genes), intent(in) :: genes
            !! Gene expression matrix (n_tissues × n_genes), column-major
        real(real64), dimension(n_tissues, n_families), intent(in) :: centroids
            !! Family centroid matrix (n_tissues × n_families), column-major, `DISTANCE_SENTINEL` for unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! M_GENE_TO_FAM_DOC(genes)
        real(real64), dimension(n_genes), intent(out) :: distances
            !! Output distances array

        integer(int32) :: i_gene, family_idx

        do concurrent (i_gene = 1:n_genes) local(family_idx) shared(gene_to_fam, n_families, distances, genes, centroids, n_tissues)
            family_idx = gene_to_fam(i_gene)
            if (family_idx < 1 .or. family_idx > n_families) then
                distances(i_gene) = DISTANCE_SENTINEL  ! Error indicator
                cycle
            end if
            call euclidean_distance_helper(genes(:, i_gene), centroids(:, family_idx), n_tissues, distances(i_gene))
        end do
    end subroutine distance_to_centroid_helper
end module tox_euclidean_distance

!> C wrapper for euclidean_distance.
!| Exposes euclidean_distance to C via iso_c_binding types.
pure subroutine euclidean_distance_c(vec1, vec2, n_elements, result, ierr) bind(C, name="euclidean_distance_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_euclidean_distance, only: euclidean_distance
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_elements
        !! Dimension of both vectors
    real(c_double), dimension(n_elements), intent(in), target :: vec1
        !! First expression vector
    real(c_double), dimension(n_elements), intent(in), target :: vec2
        !! Second expression vector
    real(c_double), intent(out), target :: result
        !! Output scalar distance
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_elements)
    M_CHECK_NON_NULL(vec1)
    M_CHECK_NON_NULL(vec2)
    M_CHECK_NON_NULL(result)

    call euclidean_distance(vec1, vec2, n_elements, result, ierr)
end subroutine euclidean_distance_c

!> C wrapper for distance_to_centroid.
!| Exposes distance_to_centroid to C via iso_c_binding types.
pure subroutine distance_to_centroid_c(n_genes, n_families, genes, centroids, &
                                       gene_to_fam, distances, n_elements, ierr) bind(C, name="distance_to_centroid_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_euclidean_distance, only: distance_to_centroid
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    integer(c_int), intent(in), target :: n_families
        !! Total number of gene families
    integer(c_int), intent(in), target :: n_elements
        !! Expression vector dimension
    real(c_double), dimension(n_elements, n_genes), intent(in), target :: genes
        !! Gene expression matrix (n_elements × n_genes), column-major
    real(c_double), dimension(n_elements, n_genes), intent(in), target :: centroids
        !! Family centroid matrix (n_elements × n_families), column-major
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! M_GENE_TO_FAM_DOC(genes)
    real(c_double), dimension(n_genes), intent(out), target :: distances
        !! Output distances array
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(n_elements)
    M_CHECK_NON_NULL(genes)
    M_CHECK_NON_NULL(centroids)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(distances)

    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                              gene_to_fam, distances, n_elements, ierr)
end subroutine distance_to_centroid_c
