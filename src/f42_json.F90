#include "macros.h"

!> This module is for JSON-specification compliant serialization of any data. In future, deserialization may be added.
module f42_json
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none

    private

    public :: serialize_json_array, serialize_json_object

    !> The core value wrapper
    !!
    !! If `value` is unassigned or of unsupported type, the result will be `null`
    type, public :: json_value
        class(*), pointer :: value => null()
            !! Supported types: `integer(int32)`, `real(real64)`, `logical`, `complex(real64)`, [[json_array(type)]], [[json_object(type)]], [[json_value(type)]]
            !!
            !! The support for [[json_value(type)]] is just for simplicity but not recommended for use, as it increases recursion depth unnecessarily.
    end type json_value

    !> Wrapper type for JSON Arrays
    !|
    !| If `elements` is unassigned, the result will be an empty array `[]`.
    type, public :: json_array
        class(*), dimension(:), pointer :: elements => null()
            !! Supported types: `integer(int32)`, `real(real64)`, `logical`, `complex(real64)`, [[json_array(type)]], [[json_object(type)]], [[json_value(type)]]
    end type json_array

    !> Wrapper type for JSON objects.
    !|
    !| - holds the key-value pairs as two separate arrays, one for `keys`, one for `values`
    !| - if the array sizes differ, the lower size is used, so only true pairs are serialized
    !| - if at least one member is not assigned to an array, the result will be an empty object `{}`
    type, public :: json_object
        character(len=:), dimension(:), pointer :: keys => null()
            !! Array of the keys of the key-value pairs
        class(*), dimension(:), pointer :: values => null()
            !! Array of the values of the key-value pairs.
            !!
            !! Supported types: `integer(int32)`, `real(real64)`, `logical`, `complex(real64)`, [[json_array(type)]], [[json_object(type)]], [[json_value(type)]]
    end type json_object
contains

    !> Serializes a real number as JSON. `Infinity` and `NaN` result in `null`
    subroutine serialize_real(real_num, unit)
        real(real64), intent(in) :: real_num
            !! real number to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to

        if (ieee_is_nan(real_num) .or. .not. ieee_is_finite(real_num)) then
            write (unit, "('null')", advance="no")
        else
            write (unit, "(ES24.16E3)", advance="no") real_num
        end if
    end subroutine serialize_real

    !> Serializes any scalar value as JSON.
    !! All intrinsics `integer(int32)`, `real(real64)`, `logical`, `complex(real64)` are supported, everything else results in `null`.
    !!
    !! `complex(real64)` will be serialized as array of the two components: [real, imag]
    subroutine serialize_scalar(scalar, unit)
        class(*), intent(in) :: scalar
            !! Scalar value to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to

        select type(scalar)
            type is (integer(int32))
                write (unit, "(I0)", advance="no") scalar
            type is (real(real64))
                call serialize_real(scalar, unit)
            type is (logical)
                write (unit, "(A)", advance="no") trim(merge("true ", "false", scalar))
            type is (character(*))
                write (unit, "('""', A, '""')", advance="no") trim(scalar)
            type is (complex(real64))
                ! Serialize as array [r, i]
                write (unit, "('[')", advance="no")
                call serialize_real(real(scalar), unit)
                write (unit, "(',')", advance="no")
                call serialize_real(aimag(scalar), unit)
                write (unit, "(']')", advance="no")
            class default
                write (unit, "('null')", advance="no")
        end select
    end subroutine serialize_scalar

    !> Serializes a [[json_array(type)]] and writes it to the passed unit
    recursive subroutine serialize_array(json_arr, unit, depth, max_depth)
        type(json_array), intent(in) :: json_arr
            !! JSON Array to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(inout) :: depth
            !! Current depth of traversion
        integer(int32), intent(in) :: max_depth
            !! maximum recursion depth for traversion of `json_arr`

        write (unit, "('[')", advance="no")
        if (associated(json_arr%elements)) then
            block
                integer(int32) :: n_elements, i_element

                n_elements = size(json_arr%elements, dim=1, kind=int32)
                do i_element = 1, n_elements
                    call serialize_json_value(json_arr%elements(i_element), unit, depth, max_depth)
                    if (i_element < n_elements) write (unit, "(',')", advance="no")
                end do
            end block
        end if
        write (unit, "(']')", advance="no")
    end subroutine serialize_array

    !> Serializes any supported type, else `null`
    recursive subroutine serialize_json_value(element, unit, depth, max_depth)
        class(*), intent(in) :: element
            !! JSON value to serialize, can be either [[json_value(type)]] or any type supported by [[json_value(type)]]
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(inout) :: depth
            !! Current depth of traversion
        integer(int32), intent(in) :: max_depth
            !! maximum recursion depth for traversion of `element`

        if (depth < max_depth) then
            depth = depth + 1
            select type(element)
                type is (json_array)
                    call serialize_array(element, unit, depth, max_depth)
                type is (json_object)
                    call serialize_object(element, unit, depth, max_depth)
                type is (json_value)
                    if (associated(element%value)) then
                        call serialize_json_value(element%value, unit, depth, max_depth)
                    else
                        write (unit, "('null')", advance="no")
                    end if
                class default
                    call serialize_scalar(element, unit)
            end select
            depth = depth - 1
        else
            write (unit, "('null')", advance="no")
        end if
    end subroutine serialize_json_value

    !> Serializes a key-value pair
    recursive subroutine serialize_key_value_pair(key, value, unit, depth, max_depth)
        character(len=*), intent(in) :: key
            !! Key of the pair
        class(*), intent(in) :: value
            !! Value of the pair
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(inout) :: depth
            !! Current depth of traversion
        integer(int32), intent(in) :: max_depth
            !! maximum recursion depth that `depth` must not exceed

        write (unit, "('""', A, '"":')", advance="no") trim(key)
        call serialize_json_value(value, unit, depth, max_depth)
    end subroutine serialize_key_value_pair

    !> Serializes a [[json_object(type)]] and writes it to the passed unit
    recursive subroutine serialize_object(json_obj, unit, depth, max_depth)
        type(json_object), intent(in) :: json_obj
            !! JSON Object to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(inout) :: depth
            !! Current depth of traversion
        integer(int32), intent(in) :: max_depth
            !! maximum recursion depth for traversion of `json_obj`

        write (unit, "('{')", advance="no")
        if (associated(json_obj%keys) .and. associated(json_obj%values)) then
            block
                integer(int32) :: n_entries, i_entry

                n_entries = min(size(json_obj%keys, dim=1, kind=int32), size(json_obj%values, dim=1, kind=int32))
                do i_entry = 1, n_entries
                    call serialize_key_value_pair(json_obj%keys(i_entry), json_obj%values(i_entry), unit, depth, max_depth)
                    if (i_entry < n_entries) write (unit, "(',')", advance="no")
                end do
            end block
        end if
        write (unit, "('}')", advance="no")
    end subroutine serialize_object

    !> Serializes a [[json_object(type)]] and writes it to the passed unit
    subroutine serialize_json_object(json_obj, unit, max_depth)
        type(json_object), intent(in) :: json_obj
            !! JSON Object to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(in), optional :: max_depth
            !! maximum recursion depth for traversion of `json_obj`, default: 20

        integer(int32) :: depth, actual_max_depth

        M_DEFAULT_VAL(max_depth, actual_max_depth, 20_int32)

        depth = 0_int32
        call serialize_object(json_obj, unit, depth, actual_max_depth)
    end subroutine serialize_json_object

    !> Serializes a [[json_array(type)]] and writes it to the passed unit
    subroutine serialize_json_array(json_arr, unit, max_depth)
        type(json_array), intent(in) :: json_arr
            !! JSON Array to serialize
        integer(int32), intent(in) :: unit
            !! unit of the file to write to
        integer(int32), intent(in), optional :: max_depth
            !! maximum recursion depth for traversion of `json_arr`, default: 20

        integer(int32) :: depth, actual_max_depth

        M_DEFAULT_VAL(max_depth, actual_max_depth, 20_int32)

        depth = 0_int32
        call serialize_array(json_arr, unit, depth, actual_max_depth)
    end subroutine serialize_json_array

    subroutine serialize_tox_data_as_flyer_json(tissues, family_ids, centroids, gene_ids, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types, ierr)
        use tox_errors, only: is_err, set_ok, set_err, ERR_INVALID_INPUT, validate_dimension_size

        character(len=*), dimension(:), contiguous, intent(in), target :: tissues
        character(len=*), dimension(:), contiguous, intent(in), target :: family_ids
        character(len=*), dimension(:), contiguous, intent(in), target :: gene_ids
        character(len=*), dimension(:), contiguous, intent(in), target :: gene_types
        character(len=*), dimension(:), contiguous, intent(in), target :: gene_species
        real(real64), dimension(:, :), contiguous, intent(in), target :: centroids
        real(real64), dimension(:, :), contiguous, intent(in), target :: genes
        integer(int32), dimension(:), contiguous, intent(in) :: gene_to_fam
        integer(int32), dimension(:), contiguous, intent(in), target :: sorted_gene_to_fam_perm
        logical, dimension(:), contiguous, intent(in), target :: gene_outliers
        integer(int32), intent(out) :: ierr

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

        integer(int32) :: i_family, i_gene, n_families, n_tissues, n_genes, first_fam_member, last_fam_member

        call set_ok(ierr)

        n_tissues = size(tissues, dim=1, kind=int32)
        n_genes = size(genes, dim=1, kind=int32)
        n_families = size(family_ids, dim=1, kind=int32)
        
        call validate_dimension_size(n_tissues, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_dimension_size(n_genes, ierr)
        if (size(gene_types, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(gene_species, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(gene_to_fam, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(sorted_gene_to_fam_perm, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(gene_outliers, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(centroids, 1, int32) /= n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (size(centroids, 2, int32) /= n_families) call set_err(ierr, ERR_INVALID_INPUT)

        if (is_err(ierr)) return

        M_ALLOCATE(character(len=8) :: flyer_json_keys(3))
        M_ALLOCATE(flyer_json_values(3))

        M_ALLOCATE(family_data(n_families))
        M_ALLOCATE(character(len=12) :: family_data_keys(3))
        M_ALLOCATE(family_data_values(size(family_data_keys, 1), n_families))
        M_ALLOCATE(family_data_gene_indices(n_families))
        M_ALLOCATE(family_data_centroids(n_families))

        M_ALLOCATE(gene_data(n_genes))
        M_ALLOCATE(character(len=11) :: gene_data_keys(6))
        M_ALLOCATE(gene_data_values(size(gene_data_keys, 1), n_genes))
        M_ALLOCATE(gene_data_coordinates(n_genes))

        flyer_json_keys = ["tissues ", "families", "genes   "]
        flyer_json_values(1)%elements => tissues
        flyer_json_values(2)%elements => family_data
        flyer_json_values(3)%elements => gene_data

        first_fam_member = 1
        do while (gene_to_fam(sorted_gene_to_fam_perm(first_fam_member)) < 1)
            first_fam_member = first_fam_member + 1
        end do

        flyer_json_keys = ["tissues ", "families", "genes   "]
        do i_family = 1, n_families
            family_data(i_family)%keys => family_data_keys
            family_data(i_family)%values => family_data_values(:, i_family)

            family_data_values(1, i_family)%value => family_ids(i_family)

            if (first_fam_member <= n_genes) then
                if (gene_to_fam(sorted_gene_to_fam_perm(first_fam_member)) == i_family) then
                    last_fam_member = first_fam_member
                    do while (gene_to_fam(sorted_gene_to_fam_perm(last_fam_member)) == i_family)
                        last_fam_member = last_fam_member + 1
                    end do
                    family_data_values(2, i_family)%value => family_data_gene_indices(i_family)
                    family_data_gene_indices(i_family)%elements => sorted_gene_to_fam_perm(first_fam_member:last_fam_member)
                    first_fam_member = last_fam_member + 1
                end if
            end if

            family_data_values(3, i_family)%value => family_data_centroids(i_family)
            family_data_centroids(i_family)%elements => centroids(:, i_family)
        end do

        gene_data_keys = ["coordinates", "id         ", "family     ", "species    ", "is_outlier ", "type       "]
        do i_gene = 1, n_genes
            gene_data_values(1, i_gene)%value => gene_data_coordinates(i_gene)
            gene_data_coordinates(i_gene)%elements => genes(:, i_gene)

            gene_data_values(2, i_gene)%value => gene_ids(i_gene)
            if (gene_to_fam(i_gene) > 0) then
                gene_data_values(3, i_gene)%value => family_ids(gene_to_fam(i_gene))
            end if
            gene_data_values(4, i_gene)%value => gene_species(i_gene)
            gene_data_values(5, i_gene)%value => gene_outliers(i_gene)
            gene_data_values(6, i_gene)%value => gene_types(i_gene)
        end do
    end subroutine serialize_tox_data_as_flyer_json
end module f42_json
