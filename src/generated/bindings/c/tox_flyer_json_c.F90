#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_flyer_json(module)]]
!| Serialization of tox analysis results into the JSON format consumed by the tox_flyer viewer.
!| Builds a [[f42_json(module)]] document model from the tox data arrays and writes it out. This
!| module is the tox-domain boundary on top of the generic [[f42_json(module)]] serializer.
module tox_flyer_json_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: serialize_tox_data_as_flyer_json_c

contains

    !> summary: C-wrapper for [[tox_flyer_json(module):serialize_tox_data_as_flyer_json(subroutine)]]
    subroutine serialize_tox_data_as_flyer_json_c(&
            filename,&
            filename_strlen,&
            tissues,&
            tissues_strlen,&
            n_tissues,&
            family_ids,&
            family_ids_strlen,&
            n_families,&
            centroids,&
            gene_ids,&
            gene_ids_strlen,&
            n_genes,&
            genes,&
            gene_to_fam,&
            sorted_gene_to_fam_perm,&
            gene_outliers,&
            gene_species,&
            gene_species_strlen,&
            gene_types,&
            gene_types_strlen,&
            ierr&
        ) bind(C, name="serialize_tox_data_as_flyer_json_c")
        use tox_flyer_json, only: serialize_tox_data_as_flyer_json

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: tissues_strlen
            !! length of the strings in `tissues`
        integer(c_int), intent(in), target :: n_tissues
            !! Number of tissues
        integer(c_int), intent(in), target :: family_ids_strlen
            !! length of the strings in `family_ids`
        integer(c_int), intent(in), target :: n_families
            !! Number of families
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: gene_species_strlen
            !! length of the strings in `gene_species`
        integer(c_int), intent(in), target :: gene_types_strlen
            !! length of the strings in `gene_types`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file to write the output to
        character(len=1, kind=c_char), dimension(tissues_strlen, n_tissues), intent(in), target :: tissues
            !! Tissue identifiers
        character(len=1, kind=c_char), dimension(family_ids_strlen, n_families), intent(in), target :: family_ids
            !! Family identifiers
        real(c_double), dimension(n_tissues, n_families), intent(in), target :: centroids
            !! Centroid data
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_genes), intent(in), target :: gene_ids
            !! Gene identifiers
        real(c_double), dimension(n_tissues, n_genes), intent(in), target :: genes
            !! Gene data
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Gene index to Family index mapping
        integer(c_int), dimension(n_genes), intent(in), target :: sorted_gene_to_fam_perm
            !! Permutation vector that sorts `gene_to_fam`
        logical(c_bool), dimension(n_genes), intent(in), target :: gene_outliers
            !! Specifies if a gene is an outlier
        character(len=1, kind=c_char), dimension(gene_species_strlen, n_genes), intent(in), target :: gene_species
            !! Species name per gene
        character(len=1, kind=c_char), dimension(gene_types_strlen, n_genes), intent(in), target :: gene_types
            !! Gene type string (ortholog/paralog)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: filename_f
        character(len=:), allocatable, dimension(:) :: tissues_f
        character(len=:), allocatable, dimension(:) :: family_ids_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        logical, dimension(n_genes) :: gene_outliers_f
        character(len=:), allocatable, dimension(:) :: gene_species_f
        character(len=:), allocatable, dimension(:) :: gene_types_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(tissues_strlen)
        M_CHECK_NON_NULL(n_tissues)
        M_CHECK_NON_NULL(family_ids_strlen)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(gene_species_strlen)
        M_CHECK_NON_NULL(gene_types_strlen)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(tissues, tissues_strlen * n_tissues)
        M_CHECK_ARRAY_NON_NULL(family_ids, family_ids_strlen * n_families)
        M_CHECK_ARRAY_NON_NULL(centroids, n_tissues * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_genes)
        M_CHECK_ARRAY_NON_NULL(genes, n_tissues * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_gene_to_fam_perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_outliers, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_species, gene_species_strlen * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_types, gene_types_strlen * n_genes)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(tissues, tissues_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(family_ids, family_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        gene_outliers_f = gene_outliers
        call c_char_2d_as_string(gene_species, gene_species_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_types, gene_types_f, ierr)
        if (is_err(ierr)) return

        call serialize_tox_data_as_flyer_json(&
            filename = filename_f,&
            tissues = tissues_f,&
            n_tissues = n_tissues,&
            family_ids = family_ids_f,&
            n_families = n_families,&
            centroids = centroids,&
            gene_ids = gene_ids_f,&
            n_genes = n_genes,&
            genes = genes,&
            gene_to_fam = gene_to_fam,&
            sorted_gene_to_fam_perm = sorted_gene_to_fam_perm,&
            gene_outliers = gene_outliers_f,&
            gene_species = gene_species_f,&
            gene_types = gene_types_f,&
            ierr = ierr&
        )
    end subroutine serialize_tox_data_as_flyer_json_c

end module tox_flyer_json_c
#endif
