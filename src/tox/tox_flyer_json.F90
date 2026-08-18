#include <src/macros.h>

!> Serialization of tox analysis results into the JSON format consumed by the tox_flyer viewer.
!| Builds a [[f42_json(module)]] document model from the tox data arrays and writes it out. This
!| module is the tox-domain boundary on top of the generic [[f42_json(module)]] serializer.
module tox_flyer_json
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_json, only: json_object, json_array, json_value, serialize_json_object
    implicit none

    private

    public :: serialize_tox_data_as_flyer_json
contains

    !> Serializes tox related data to JSON, compatible with the tox_flyer
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


subroutine serialize_tox_data_as_flyer_json_c(filename, filename_len, tissues, tissue_len, n_tissues, family_ids, family_id_len, n_families, centroids, gene_ids, gene_id_len, n_genes, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_species_len, gene_types, gene_type_len, ierr) bind(C, name="serialize_tox_data_as_flyer_json_c")

    use tox_flyer_json, only: serialize_tox_data_as_flyer_json
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_double
    use tox_conversions, only: c_char_2d_as_string, c_char_1d_as_string, c_int_as_logical
    use tox_errors, only: set_ok, validate_dimension_size, is_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_tissues
        !! Number of tissues
    integer(c_int), intent(in), target :: n_families
        !! Number of families
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: filename_len
        !! String Length of `filename`
    integer(c_int), intent(in), target :: tissue_len
        !! String Length of `tissues`
    integer(c_int), intent(in), target :: family_id_len
        !! String Length of `family_ids`
    integer(c_int), intent(in), target :: gene_id_len
        !! String Length of `gene_ids`
    integer(c_int), intent(in), target :: gene_type_len
        !! String Length of `gene_types`
    integer(c_int), intent(in), target :: gene_species_len
        !! String Length of `gene_species`
    character(kind=c_char), dimension(filename_len), intent(in), target :: filename
        !! Name of the file to write the output to
    character(kind=c_char), dimension(tissue_len, n_tissues), intent(in), target :: tissues
        !! Tissue identifiers
    character(kind=c_char), dimension(family_id_len, n_families), intent(in), target :: family_ids
        !! Family identifiers
    character(kind=c_char), dimension(gene_id_len, n_genes), intent(in), target :: gene_ids
        !! Gene identifiers
    character(kind=c_char), dimension(gene_type_len, n_genes), intent(in), target :: gene_types
        !! Gene type string (ortholog/paralog)
    character(kind=c_char), dimension(gene_species_len, n_genes), intent(in), target :: gene_species
        !! Species name per gene
    real(c_double), dimension(n_tissues, n_families), intent(in), target :: centroids
        !! Centroid data
    real(c_double), dimension(n_tissues, n_genes), intent(in), target :: genes
        !! Gene data
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene index to Family index mapping
    integer(c_int), dimension(n_genes), intent(in), target :: sorted_gene_to_fam_perm
        !! Permutation vector that sorts `gene_to_fam`
    integer(c_int), dimension(n_genes), intent(in), target :: gene_outliers
        !! Specifies if a gene is an outlier
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: f_filename
    character(len=:), allocatable :: f_tissues(:)
    character(len=:), allocatable :: f_family_ids(:)
    character(len=:), allocatable :: f_gene_ids(:)
    character(len=:), allocatable :: f_gene_species(:)
    character(len=:), allocatable :: f_gene_types(:)
    logical, allocatable :: f_gene_outliers(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_tissues)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(filename_len)
    M_CHECK_NON_NULL(tissue_len)
    M_CHECK_NON_NULL(family_id_len)
    M_CHECK_NON_NULL(gene_id_len)
    M_CHECK_NON_NULL(gene_type_len)
    M_CHECK_NON_NULL(gene_species_len)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(tissues)
    M_CHECK_NON_NULL(family_ids)
    M_CHECK_NON_NULL(gene_ids)
    M_CHECK_NON_NULL(gene_types)
    M_CHECK_NON_NULL(gene_species)
    M_CHECK_NON_NULL(centroids)
    M_CHECK_NON_NULL(genes)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(sorted_gene_to_fam_perm)
    M_CHECK_NON_NULL(gene_outliers)

    call c_char_1d_as_string(filename, f_filename, ierr)
    if(is_err(ierr)) return
    call c_char_2d_as_string(tissues, f_tissues, ierr)
    if(is_err(ierr)) return
    call c_char_2d_as_string(family_ids, f_family_ids, ierr)
    if(is_err(ierr)) return
    call c_char_2d_as_string(gene_ids, f_gene_ids, ierr)
    if(is_err(ierr)) return
    call c_char_2d_as_string(gene_species, f_gene_species, ierr)
    if(is_err(ierr)) return
    call c_char_2d_as_string(gene_types, f_gene_types, ierr)
    if(is_err(ierr)) return

    M_ALLOCATE(f_gene_outliers(n_genes))
    call c_int_as_logical(gene_outliers, f_gene_outliers)

    call serialize_tox_data_as_flyer_json( &
        f_filename, &
        f_tissues, n_tissues, &
        f_family_ids, n_families, &
        centroids, &
        f_gene_ids, n_genes, &
        genes, &
        gene_to_fam, &
        sorted_gene_to_fam_perm, &
        f_gene_outliers, &
        f_gene_species, &
        f_gene_types, &
        ierr )

end subroutine serialize_tox_data_as_flyer_json_c
