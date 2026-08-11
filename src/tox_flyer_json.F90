#include <src/macros.h>

!> Serialization of tox analysis results into the JSON format consumed by the tox_flyer viewer.
!| Builds a [[f42_json(module)]] document model from the tox data arrays and writes it out. This
!| module is the tox-domain boundary on top of the generic [[f42_json(module)]] serializer.
module tox_flyer_json
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_json, only: json_object, json_array, json_value, serialize_json_object
    M_IMPLICIT_NONE

    private

    public :: serialize_tox_data_as_flyer_json
contains

    !> M_EXPORT_C
    !| summary: Serializes tox related data to JSON, compatible with the tox_flyer
    !| AUTHOR_FRANZ_ERIC_SILL
    subroutine serialize_tox_data_as_flyer_json(filename, tissues, n_tissues, family_ids, n_families, centroids, gene_ids, n_genes, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types, ierr)
        use tox_errors, only: is_err, set_ok, set_err, ERR_INVALID_INPUT, ERR_ALLOC_FAIL, validate_dimension_size, validate_all_in_range_int
        use f42_utils, only: open_file

        integer(int32), intent(in) :: n_tissues
            !! Number of tissues
        integer(int32), intent(in) :: n_families
            !! Number of families
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        character(len=*), intent(in) :: filename
            !! Name of the file to write the output to
        character(len=*), dimension(n_tissues), intent(in), target :: tissues
            !! Tissue identifiers
        character(len=*), dimension(n_families), intent(in), target :: family_ids
            !! Family identifiers
        character(len=*), dimension(n_genes), intent(in), target :: gene_ids
            !! Gene identifiers
        character(len=*), dimension(n_genes), intent(in), target :: gene_types
            !! Gene type string (ortholog/paralog)
        character(len=*), dimension(n_genes), intent(in), target :: gene_species
            !! Species name per gene
        real(real64), dimension(n_tissues, n_families), intent(in), target :: centroids
            !! Centroid data
        real(real64), dimension(n_tissues, n_genes), intent(in), target :: genes
            !! Gene data
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene index to Family index mapping
        integer(int32), dimension(n_genes), intent(in), target :: sorted_gene_to_fam_perm
            !! Permutation vector that sorts `gene_to_fam`
        logical, dimension(n_genes), intent(in), target :: gene_outliers
            !! Specifies if a gene is an outlier
        integer(int32), intent(out) :: ierr
            !! Error code

        type(json_object) :: flyer_json
        character(len=:), dimension(:), allocatable, target :: flyer_json_keys
        type(json_array), dimension(:), allocatable, target :: flyer_json_values

        type(json_object), dimension(:), allocatable, target :: family_data
        character(len=:), dimension(:), allocatable, target :: family_data_keys
        type(json_value), dimension(:, :), allocatable, target :: family_data_values
        type(json_array), dimension(:), allocatable, target :: family_data_gene_indices
        type(json_array), dimension(:), allocatable, target :: family_data_centroids

        type(json_object), dimension(:), allocatable, target :: gene_data
        character(len=:), dimension(:), allocatable, target :: gene_data_keys
        type(json_value), dimension(:, :), allocatable, target :: gene_data_values
        type(json_array), dimension(:), allocatable, target :: gene_data_coordinates

        integer(int32) :: i_family, i_gene, first_fam_member, last_fam_member, unit

        call set_ok(ierr)

        call validate_dimension_size(n_tissues, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_all_in_range_int(sorted_gene_to_fam_perm, n_genes, ierr, min=1_int32, max=n_genes)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, min=0_int32, max=n_families)

        if (is_err(ierr)) return

        ! create gene json array [{coordinates:..., id:..., is_outlier:..., ...}, ...]
        M_ALLOCATE(gene_data(n_genes))
        M_ALLOCATE(character(len=11) :: gene_data_keys(6))
        M_ALLOCATE(gene_data_values(size(gene_data_keys, 1), n_genes))
        M_ALLOCATE(gene_data_coordinates(n_genes))

        ! create root json object {tissues:..., families:..., genes:...}
        M_ALLOCATE(character(len=8) :: flyer_json_keys(3))
        M_ALLOCATE(flyer_json_values(3))
        flyer_json_keys = ["tissues ", "families", "genes   "]
        flyer_json_values(1)%elements => tissues
        flyer_json_values(3)%elements => gene_data
        flyer_json%keys => flyer_json_keys
        flyer_json%values => flyer_json_values

        ! create family json array [{family:..., gene_indices:..., centroid:..., ...}, ...]
        M_ALLOCATE(family_data(n_families))
        flyer_json_values(2)%elements => family_data
        M_ALLOCATE(character(len=12) :: family_data_keys(3))

        M_ALLOCATE(family_data_values(size(family_data_keys, 1), n_families))
        M_ALLOCATE(family_data_gene_indices(n_families))
        M_ALLOCATE(family_data_centroids(n_families))

        ! find index of first family (skip unassigned genes)
        first_fam_member = 1
        do i_gene = 1, n_genes
            if (gene_to_fam(sorted_gene_to_fam_perm(first_fam_member)) > 0) exit
            first_fam_member = first_fam_member + 1
        end do

        ! create the family object per family
        family_data_keys = ["family      ", "gene_indices", "centroid    "]
        do i_family = 1, n_families

            ! point to family id
            family_data_values(1, i_family)%value => family_ids(i_family)

            ! point to gene indices array
            family_data_values(2, i_family)%value => family_data_gene_indices(i_family)

            ! find last gene of the current family and point to the slice of gene indices
            if (first_fam_member <= n_genes) then
                last_fam_member = first_fam_member
                do while (gene_to_fam(sorted_gene_to_fam_perm(last_fam_member)) == i_family)
                    last_fam_member = last_fam_member + 1
                    if (last_fam_member > n_genes) exit
                end do
                last_fam_member = last_fam_member - 1
                ! point to gene indices
                family_data_gene_indices(i_family)%elements => sorted_gene_to_fam_perm(first_fam_member:last_fam_member)
                first_fam_member = last_fam_member + 1
            end if

            ! point to centroid
            family_data_values(3, i_family)%value => family_data_centroids(i_family)
            family_data_centroids(i_family)%elements => centroids(:, i_family)


            ! assign to family_data array
            family_data(i_family)%keys => family_data_keys
            family_data(i_family)%values => family_data_values(:, i_family)
        end do

        ! 3. create array of gene objects
        gene_data_keys = ["coordinates", "id         ", "family     ", "species    ", "is_outlier ", "type       "]
        do i_gene = 1, n_genes

            ! point to coordinates
            gene_data_values(1, i_gene)%value => gene_data_coordinates(i_gene)
            gene_data_coordinates(i_gene)%elements => genes(:, i_gene)

            ! point to gene id
            gene_data_values(2, i_gene)%value => gene_ids(i_gene)

            ! point to family id
            if (gene_to_fam(i_gene) > 0) then
                gene_data_values(3, i_gene)%value => family_ids(gene_to_fam(i_gene))
            end if

            ! point to species name
            gene_data_values(4, i_gene)%value => gene_species(i_gene)

            ! point to logical outlier specifier
            gene_data_values(5, i_gene)%value => gene_outliers(i_gene)

            ! point to gene type specifier
            gene_data_values(6, i_gene)%value => gene_types(i_gene)


            ! assign to gene_data array
            gene_data(i_gene)%keys => gene_data_keys
            gene_data(i_gene)%values => gene_data_values(:, i_gene)
        end do

        ! open file for formatted data
        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return
        call serialize_json_object(flyer_json, unit)
        close(unit)
    end subroutine serialize_tox_data_as_flyer_json

end module tox_flyer_json
