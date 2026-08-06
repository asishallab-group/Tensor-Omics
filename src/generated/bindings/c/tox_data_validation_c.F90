#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_validation(module)]]
!| Semantic validation of TensorOmics data sets (dimensions, ID uniqueness, value ranges,
!| and cross-array consistency).
!|
!| Complements [[tox_errors(module)]]'s generic per-argument validators
!| (`validate_dimension_size`, `validate_in_range_int/real`, ...) with checks specific to the
!| TensorOmics gene/family/expression/centroid/shift-vector data model, such as verifying that
!| `shift_vectors` was actually derived from `expression_vectors` and `family_centroids` as
!| expected. [[tox_data_validation(module):validate_all_data(subroutine)]] runs the full suite.
module tox_data_validation_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_2d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: validate_data_structure_c
    public :: validate_gene_to_family_mapping_c
    public :: validate_expression_data_c
    public :: validate_family_centroids_c
    public :: validate_shift_vectors_c
    public :: validate_string_array_uniqueness_c
    public :: validate_all_data_c

contains

    !> summary: C-wrapper for [[tox_data_validation(module):validate_data_structure(subroutine)]]
    subroutine validate_data_structure_c(&
            n_genes,&
            n_families,&
            n_samples,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            gene_family_ids,&
            gene_family_ids_strlen,&
            n_gene_family_ids_elements,&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            expression_vectors,&
            n_expression_vectors_elements_dim_1,&
            n_expression_vectors_elements_dim_2,&
            family_centroids,&
            n_family_centroids_elements_dim_1,&
            n_family_centroids_elements_dim_2,&
            shift_vectors,&
            n_shift_vectors_elements_dim_1,&
            n_shift_vectors_elements_dim_2,&
            ierr&
        ) bind(C, name="validate_data_structure_c")
        use tox_data_validation, only: validate_data_structure

        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        integer(c_int), intent(in), target :: gene_family_ids_strlen
            !! length of the strings in `gene_family_ids`
        integer(c_int), intent(in), target :: n_gene_family_ids_elements
            !! number of elements in `gene_family_ids`
        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !! 1. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !! 2. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !! 1. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !! 2. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_1
            !! 1. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_2
            !! 2. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_genes
            !! Expected number of genes
        integer(c_int), intent(in), target :: n_families
            !! Expected number of families
        integer(c_int), intent(in), target :: n_samples
            !! Expected number of samples
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(in), target :: gene_ids
            !! Gene ids
        character(len=1, kind=c_char), dimension(gene_family_ids_strlen, n_gene_family_ids_elements), intent(in), target :: gene_family_ids
            !! Gene family ids
        integer(c_int), dimension(n_gene_to_fam_elements), intent(in), target :: gene_to_fam
            !! gene to family mapping
        real(c_double), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), intent(in), target :: expression_vectors
            !! Expression vectors
        real(c_double), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), intent(in), target :: family_centroids
            !! Family centroids
        real(c_double), dimension(n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2), intent(in), target :: shift_vectors
            !! Shift vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable, dimension(:) :: gene_family_ids_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(gene_family_ids_strlen)
        M_CHECK_NON_NULL(n_gene_family_ids_elements)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_gene_ids_elements)
        M_CHECK_ARRAY_NON_NULL(gene_family_ids, gene_family_ids_strlen * n_gene_family_ids_elements)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_expression_vectors_elements_dim_1 * n_expression_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_family_centroids_elements_dim_1 * n_family_centroids_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(shift_vectors, n_shift_vectors_elements_dim_1 * n_shift_vectors_elements_dim_2)

        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_family_ids, gene_family_ids_f, ierr)
        if (is_err(ierr)) return

        call validate_data_structure(&
            n_genes = n_genes,&
            n_families = n_families,&
            n_samples = n_samples,&
            gene_ids = gene_ids_f,&
            gene_family_ids = gene_family_ids_f,&
            gene_to_fam = gene_to_fam,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            shift_vectors = shift_vectors,&
            ierr = ierr&
        )
    end subroutine validate_data_structure_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_gene_to_family_mapping(subroutine)]]
    subroutine validate_gene_to_family_mapping_c(&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            n_families,&
            ierr&
        ) bind(C, name="validate_gene_to_family_mapping_c")
        use tox_data_validation, only: validate_gene_to_family_mapping

        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        integer(c_int), dimension(n_gene_to_fam_elements), intent(in), target :: gene_to_fam
            !! gene to family mapping
        integer(c_int), intent(in), target :: n_families
            !! number of families
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)

        call validate_gene_to_family_mapping(&
            gene_to_fam = gene_to_fam,&
            n_families = n_families,&
            ierr = ierr&
        )
    end subroutine validate_gene_to_family_mapping_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_expression_data(subroutine)]]
    subroutine validate_expression_data_c(&
            expression_vectors,&
            n_expression_vectors_elements_dim_1,&
            n_expression_vectors_elements_dim_2,&
            check_non_negative,&
            ierr&
        ) bind(C, name="validate_expression_data_c")
        use tox_data_validation, only: validate_expression_data

        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !! 1. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !! 2. dimension of `expression_vectors`
        real(c_double), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), intent(in), target :: expression_vectors
            !! Expression vectors
        logical(c_bool), intent(in), target :: check_non_negative
            !! Defines if non negative should be checked
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: check_non_negative_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(check_non_negative)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_expression_vectors_elements_dim_1 * n_expression_vectors_elements_dim_2)

        check_non_negative_f = check_non_negative

        call validate_expression_data(&
            expression_vectors = expression_vectors,&
            check_non_negative = check_non_negative_f,&
            ierr = ierr&
        )
    end subroutine validate_expression_data_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_family_centroids(subroutine)]]
    subroutine validate_family_centroids_c(&
            family_centroids,&
            n_family_centroids_elements_dim_1,&
            n_family_centroids_elements_dim_2,&
            ierr&
        ) bind(C, name="validate_family_centroids_c")
        use tox_data_validation, only: validate_family_centroids

        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !! 1. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !! 2. dimension of `family_centroids`
        real(c_double), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), intent(in), target :: family_centroids
            !! Family centroids array
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_family_centroids_elements_dim_1 * n_family_centroids_elements_dim_2)

        call validate_family_centroids(&
            family_centroids = family_centroids,&
            ierr = ierr&
        )
    end subroutine validate_family_centroids_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_shift_vectors(subroutine)]]
    subroutine validate_shift_vectors_c(&
            shift_vectors,&
            n_shift_vectors_elements_dim_1,&
            n_shift_vectors_elements_dim_2,&
            expression_vectors,&
            n_expression_vectors_elements_dim_1,&
            n_expression_vectors_elements_dim_2,&
            family_centroids,&
            n_family_centroids_elements_dim_1,&
            n_family_centroids_elements_dim_2,&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            n_samples,&
            ierr&
        ) bind(C, name="validate_shift_vectors_c")
        use tox_data_validation, only: validate_shift_vectors

        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_1
            !! 1. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_2
            !! 2. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !! 1. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !! 2. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !! 1. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !! 2. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        real(c_double), dimension(n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2), intent(in), target :: shift_vectors
            !! shift vectors
        real(c_double), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), intent(in), target :: expression_vectors
            !! expression vectors
        real(c_double), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), intent(in), target :: family_centroids
            !! family centroids
        integer(c_int), dimension(n_gene_to_fam_elements), intent(in), target :: gene_to_fam
            !! gene to family mapping
        integer(c_int), intent(in), target :: n_samples
            !! Number of samples
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_ARRAY_NON_NULL(shift_vectors, n_shift_vectors_elements_dim_1 * n_shift_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_expression_vectors_elements_dim_1 * n_expression_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_family_centroids_elements_dim_1 * n_family_centroids_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)

        call validate_shift_vectors(&
            shift_vectors = shift_vectors,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            gene_to_fam = gene_to_fam,&
            n_samples = n_samples,&
            ierr = ierr&
        )
    end subroutine validate_shift_vectors_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_string_array_uniqueness(subroutine)]]
    subroutine validate_string_array_uniqueness_c(&
            str_arr,&
            str_arr_strlen,&
            n_str_arr_elements,&
            ierr&
        ) bind(C, name="validate_string_array_uniqueness_c")
        use tox_data_validation, only: validate_string_array_uniqueness

        integer(c_int), intent(in), target :: str_arr_strlen
            !! length of the strings in `str_arr`
        integer(c_int), intent(in), target :: n_str_arr_elements
            !! number of elements in `str_arr`
        character(len=1, kind=c_char), dimension(str_arr_strlen, n_str_arr_elements), intent(in), target :: str_arr
            !! string array
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable, dimension(:) :: str_arr_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(str_arr_strlen)
        M_CHECK_NON_NULL(n_str_arr_elements)
        M_CHECK_ARRAY_NON_NULL(str_arr, str_arr_strlen * n_str_arr_elements)

        call c_char_2d_as_string(str_arr, str_arr_f, ierr)
        if (is_err(ierr)) return

        call validate_string_array_uniqueness(&
            str_arr = str_arr_f,&
            ierr = ierr&
        )
    end subroutine validate_string_array_uniqueness_c

    !> summary: C-wrapper for [[tox_data_validation(module):validate_all_data(subroutine)]]
    subroutine validate_all_data_c(&
            n_genes,&
            n_families,&
            n_samples,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            gene_family_ids,&
            gene_family_ids_strlen,&
            n_gene_family_ids_elements,&
            gene_to_fam,&
            n_gene_to_fam_elements,&
            expression_vectors,&
            n_expression_vectors_elements_dim_1,&
            n_expression_vectors_elements_dim_2,&
            family_centroids,&
            n_family_centroids_elements_dim_1,&
            n_family_centroids_elements_dim_2,&
            shift_vectors,&
            n_shift_vectors_elements_dim_1,&
            n_shift_vectors_elements_dim_2,&
            ierr,&
            check_uniqueness,&
            check_shift_consistency&
        ) bind(C, name="validate_all_data_c")
        use tox_data_validation, only: validate_all_data

        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        integer(c_int), intent(in), target :: gene_family_ids_strlen
            !! length of the strings in `gene_family_ids`
        integer(c_int), intent(in), target :: n_gene_family_ids_elements
            !! number of elements in `gene_family_ids`
        integer(c_int), intent(in), target :: n_gene_to_fam_elements
            !! number of elements in `gene_to_fam`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_1
            !! 1. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_expression_vectors_elements_dim_2
            !! 2. dimension of `expression_vectors`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !! 1. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !! 2. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_1
            !! 1. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_2
            !! 2. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of families
        integer(c_int), intent(in), target :: n_samples
            !! Number of samples
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(in), target :: gene_ids
            !! Gene ids array
        character(len=1, kind=c_char), dimension(gene_family_ids_strlen, n_gene_family_ids_elements), intent(in), target :: gene_family_ids
            !! gene family ids
        integer(c_int), dimension(n_gene_to_fam_elements), intent(in), target :: gene_to_fam
            !! gene to family mapping
        real(c_double), dimension(n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2), intent(in), target :: expression_vectors
            !! Expression vectors
        real(c_double), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), intent(in), target :: family_centroids
            !! family centroids
        real(c_double), dimension(n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2), intent(in), target :: shift_vectors
            !! shift vectors
        integer(c_int), intent(out), target :: ierr
            !! error code
        logical(c_bool), intent(in), target :: check_uniqueness
            !! Check ID arrays for uniqueness.
            !! The default value is `.true.`.
        logical(c_bool), intent(in), target :: check_shift_consistency
            !! Check consitency of shift array.
            !! The default value is `.true.`.
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable, dimension(:) :: gene_family_ids_f
        logical :: check_uniqueness_f
        logical :: check_shift_consistency_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(gene_family_ids_strlen)
        M_CHECK_NON_NULL(n_gene_family_ids_elements)
        M_CHECK_NON_NULL(n_gene_to_fam_elements)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_vectors_elements_dim_2)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_2)
        M_CHECK_NON_NULL(check_uniqueness)
        M_CHECK_NON_NULL(check_shift_consistency)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_ids_strlen * n_gene_ids_elements)
        M_CHECK_ARRAY_NON_NULL(gene_family_ids, gene_family_ids_strlen * n_gene_family_ids_elements)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_gene_to_fam_elements)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_expression_vectors_elements_dim_1 * n_expression_vectors_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_family_centroids_elements_dim_1 * n_family_centroids_elements_dim_2)
        M_CHECK_ARRAY_NON_NULL(shift_vectors, n_shift_vectors_elements_dim_1 * n_shift_vectors_elements_dim_2)

        call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(gene_family_ids, gene_family_ids_f, ierr)
        if (is_err(ierr)) return
        check_uniqueness_f = check_uniqueness
        check_shift_consistency_f = check_shift_consistency

        call validate_all_data(&
            n_genes = n_genes,&
            n_families = n_families,&
            n_samples = n_samples,&
            gene_ids = gene_ids_f,&
            gene_family_ids = gene_family_ids_f,&
            gene_to_fam = gene_to_fam,&
            expression_vectors = expression_vectors,&
            family_centroids = family_centroids,&
            shift_vectors = shift_vectors,&
            ierr = ierr,&
            check_uniqueness = check_uniqueness_f,&
            check_shift_consistency = check_shift_consistency_f&
        )
    end subroutine validate_all_data_c

end module tox_data_validation_c
#endif
