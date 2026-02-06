module tox_data_validation
    use safeguard
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
    use config, only: DEBUG
    implicit none
    private

    ! Public procedures
    public :: validate_data_structure
    public :: validate_gene_to_family_mapping
    public :: validate_expression_data
    public :: validate_family_centroids
    public :: validate_shift_vectors
    public :: check_for_nan_inf
    public :: validate_string_array_uniqueness
    public :: validate_empty_strings
    public :: validate_all_data
    
    ! Parameters for validation tolerances
    real(real64), parameter :: FLOAT_TOLERANCE = 1.0e-10_real64
    
contains

    !> Validate full data structure
    subroutine validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                                     gene_to_fam, expression_vectors, family_centroids, &
                                     shift_vectors, ierr)
        integer(int32), intent(in) :: n_genes
            !! Expected number of genes
        integer(int32), intent(in) :: n_families
            !! Expected number of families
        integer(int32), intent(in) :: n_samples
            !! Expected number of samples
        character(len=*), intent(in) :: gene_ids(:)
            !! Gene ids
        character(len=*), intent(in) :: gene_family_ids(:)
            !! Gene family ids
        integer(int32), intent(in) :: gene_to_fam(:)
            !! gene to family mapping
        real(real64), intent(in) :: expression_vectors(:,:)
            !! Expression vectors
        real(real64), intent(in) :: family_centroids(:,:)
            !! Family centroids
        real(real64), intent(in) :: shift_vectors(:,:)
            !! Shift vectors
        integer(int32), intent(out) :: ierr
            !! Error code
        
        call set_ok(ierr)
        
        ! Check basic dimensions
        if (n_genes <= 0 .or. n_families < 0 .or. n_samples <= 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Check gene_ids array
        if (size(gene_ids) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            if(DEBUG) write(*,*) 'Error: gene_ids size mismatch. Expected:', n_genes, ' Actual:', size(gene_ids)
            return
        end if
        
        ! Check gene_family_ids array
        if (size(gene_family_ids) /= n_families) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            if(DEBUG) write(*,*) 'Error: gene_family_ids size mismatch. Expected:', n_families, ' Actual:', size(gene_family_ids)
            return
        end if
        
        ! Check gene_to_fam array
        if (size(gene_to_fam) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            if(DEBUG) write(*,*) 'Error: gene_to_fam size mismatch. Expected:', n_genes, ' Actual:', size(gene_to_fam)
            return
        end if
        
        ! Check expression_vectors array
        if (size(expression_vectors, 1) /= n_samples .or. size(expression_vectors, 2) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            if(DEBUG) write(*,*) 'Error: expression_vectors size mismatch. Expected: (', n_samples, ',', n_genes, &
                       ') Actual: (', size(expression_vectors, 1), ',', size(expression_vectors, 2), ')'
            return
        end if
        
        ! Check family_centroids array if families exist
        if (n_families > 0) then
            if (size(family_centroids, 1) /= n_samples .or. size(family_centroids, 2) /= n_families) then
                call set_err_once(ierr, ERR_SIZE_MISMATCH)
                if(DEBUG) write(*,*) 'Error: family_centroids size mismatch. Expected: (', n_samples, ',', n_families, &
                           ') Actual: (', size(family_centroids, 1), ',', size(family_centroids, 2), ')'
                return
            end if
        else 
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Check shift_vectors array
        if (size(shift_vectors, 1) /= 2*n_samples .or. size(shift_vectors, 2) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            if(DEBUG) write(*,*) 'Error: shift_vectors size mismatch. Expected: (', 2*n_samples, ',', n_genes, &
                       ') Actual: (', size(shift_vectors, 1), ',', size(shift_vectors, 2), ')'
            return
        end if
        
        ! Check for empty strings in gene_ids
        call validate_empty_strings(gene_ids, "gene_ids", ierr)
        if (.not. is_ok(ierr)) return
        
        ! Check for empty strings in gene_family_ids if families exist
        if (n_families > 0) then
            call validate_empty_strings(gene_family_ids, "gene_family_ids", ierr)
            if (.not. is_ok(ierr)) return
        end if
        
    end subroutine validate_data_structure

    !> Validate gene to family mapping
    subroutine validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        integer(int32), intent(in) :: gene_to_fam(:)
            !! gene to family mapping
        integer(int32), intent(in) :: n_families
            !! number of families
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer :: i, invalid_count
        
        call set_ok(ierr)
        invalid_count = 0
        
        if (n_families == 0) then
            ! No families defined, mapping should be all zeros
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) /= 0) then
                    invalid_count = invalid_count + 1
                end if
            end do
            if (invalid_count > 0) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if(DEBUG) write(*,*) 'Error: gene_to_fam should be all zeros when no families are defined'
                RETURN
            end if
        else 
            ! Check that all gene_to_fam values are valid family indices
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) < 0 .or. gene_to_fam(i) > n_families) then
                    invalid_count = invalid_count + 1
                    if(invalid_count > 0) then
                        if(DEBUG) write(*,*) 'Error: gene_to_fam(', i, ') = ', gene_to_fam(i), &
                                   ' but valid range is 0 to ', n_families
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        return
                    end if
                end if
            end do
        end if
        
    end subroutine validate_gene_to_family_mapping

    !> Validate expresssion data
    subroutine validate_expression_data(expression_vectors, check_non_negative, ierr)
        real(real64), intent(in) :: expression_vectors(:,:)
            !! Expression vectors
        logical, intent(in) :: check_non_negative
            !! Defines if non negative should be checked
        integer(int32), intent(out) :: ierr
            !! Error code
        
        call set_ok(ierr)
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(expression_vectors, ierr)
        if (.not. is_ok(ierr)) return
        
        ! Check for negative values if requested
        if (check_non_negative) then
            if (any(expression_vectors < 0.0_real64)) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if(DEBUG) write(*,*) 'Error: Negative values found in expression data'
                return
            end if
        end if
        
    end subroutine validate_expression_data

    !> Validate the family centroids
    subroutine validate_family_centroids(family_centroids, ierr)
        real(real64), intent(in) :: family_centroids(:,:)
            !! Family centroids array
        integer(int32), intent(out) :: ierr
            !! Error code
        
        call set_ok(ierr)
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(family_centroids, ierr)
        if (.not. is_ok(ierr)) return

    end subroutine validate_family_centroids

    !> Validates shift vectors
    subroutine validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                        gene_to_fam, n_samples, ierr)
        real(real64), intent(in) :: shift_vectors(:,:)
            !! shift vectors
        real(real64), intent(in) :: expression_vectors(:,:)
            !! expression vectors
        real(real64), intent(in) :: family_centroids(:,:)
            !! family centroids
        integer(int32), intent(in) :: gene_to_fam(:)
            !! gene to family mapping
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64) :: expected_shift

        
        integer(int32) :: i, j, fam_idx, n_genes, error_count
        
        call set_ok(ierr)
        n_genes = size(expression_vectors, 2)
        error_count = 0
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(shift_vectors, ierr)
        if (.not. is_ok(ierr)) return
        
        ! Verify shift vectors structure: first d rows should be centroids, next d rows should be shifts
        do i = 1, n_genes
            fam_idx = gene_to_fam(i)
            
            if (fam_idx == 0) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                RETURN
            end if
            
            ! Check that centroid part (first d rows) matches the family centroid
            do j = 1, n_samples
                if (abs(shift_vectors(j, i) - family_centroids(j, fam_idx)) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        if(DEBUG) write(*,*) 'Error: Centroid mismatch for gene ', i, &
                                ' dimension ', j, ' expected ', family_centroids(j, fam_idx), &
                                ' got ', shift_vectors(j, i)
                    end if
                end if
            end do
            if(.not. is_ok(ierr)) return
            
            ! Check that shift part (rows d+1 to 2d) matches expression - centroid
            do j = 1, n_samples
                ! The expected shift is simply the difference between the expression vector and the centroid.
                ! This aligns with the compute_shift_vector_field subroutine.
                expected_shift = expression_vectors(j, i) - family_centroids(j, fam_idx)
                
                if (abs(shift_vectors(n_samples+j, i) - expected_shift) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        if(DEBUG) write(*,*) 'Error: Shift mismatch for gene ', i, &
                                ' dimension ', j, ' expected ', expected_shift, &
                                ' got ', shift_vectors(n_samples+j, i)
                    end if
                end if
            end do
            if(.not. is_ok(ierr)) return
        end do
        
    end subroutine validate_shift_vectors

    !> Check if an array of type real contains NaN or Inf values
    subroutine check_for_nan_inf(array, ierr)
        real(real64), intent(in) :: array(:,:)
            !! Input array
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(int32) :: i, j, nan_inf_count
        
        call set_ok(ierr)
        nan_inf_count = 0
        
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                if (is_nan(array(i, j)) .or. is_inf(array(i, j))) then
                    nan_inf_count = nan_inf_count + 1
                    if (nan_inf_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        if(DEBUG) write(*,*) 'Error: NaN/Inf found at position (', i, ',', j, ')'
                    end if
                end if
            end do
        end do
        if(.not. is_ok(ierr)) return
        
    contains

        !Possibly replace by existing asserts

        logical function is_nan(x)
            real(real64), intent(in) :: x
            is_nan = (x /= x)  ! NaN is the only value not equal to itself
        end function is_nan
        
        logical function is_inf(x)
            real(real64), intent(in) :: x
            is_inf = (x > huge(x) .or. x < -huge(x))
        end function is_inf
        
    end subroutine check_for_nan_inf

    !> Validate that no string appears more than once
    subroutine validate_string_array_uniqueness(str_arr, ierr)
        use f42_xxh3_hashmap, only: hashset_type, hashset_create, hashset_put, hashset_destroy
        character(len=*), intent(in) :: str_arr(:)
            !! string array
        integer(int32), intent(out) :: ierr
            !! Error code 
        integer(int32) :: i
        type(hashset_type) :: temp_set

        call set_ok(ierr)   
        call hashset_create(temp_set, initial_size=size(str_arr))     
        
        do i = 1, size(str_arr)
            call hashset_put(temp_set, str_arr(i), ierr)
            if (.not. is_ok(ierr)) then
                if(DEBUG) write(*,*) 'Error: Duplicate string found: ', str_arr(i)
                call hashset_destroy(temp_set)
                return
            end if
        end do
        call hashset_destroy(temp_set)
        
    end subroutine validate_string_array_uniqueness

    !> Validate that strings in an array are not empty
    subroutine validate_empty_strings(string_array, array_name, ierr)
        character(len=*), intent(in) :: string_array(:)
            !! Array of strings
        character(len=*), intent(in) :: array_name
            !! Name of the array
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(int32) :: i, empty_count
        
        call set_ok(ierr)
        empty_count = 0
        
        do i = 1, size(string_array)
            if (len_trim(string_array(i)) == 0) then
                empty_count = empty_count + 1
                if (empty_count <= 10) then
                    call set_err_once(ierr, ERR_INVALID_INPUT)
                    if(DEBUG) write(*,*) 'Error: Empty string found in ', array_name, ' at index ', i
                end if
            end if
        end do
        
    end subroutine validate_empty_strings

    !> Comprehensive validation routine, combining all checks
    subroutine validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                               gene_to_fam, expression_vectors, family_centroids, &
                               shift_vectors, ierr, check_uniqueness, check_shift_consistency)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of families
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        character(len=*), intent(in) :: gene_ids(:)
            !! Gene ids array
        character(len=*), intent(in) :: gene_family_ids(:)
            !! gene family ids
        integer(int32), intent(in) :: gene_to_fam(:)
            !! gene to family mapping
        real(real64), intent(in) :: expression_vectors(:,:)
            !! Expression vectors
        real(real64), intent(in) :: family_centroids(:,:)
            !! family centroids
        real(real64), intent(in) :: shift_vectors(:,:)
            !! shift vectors
        integer(int32), intent(out) :: ierr
            !! error code
        logical, intent(in), optional :: check_uniqueness
            !! Check ID arrays for uniqueness
        logical, intent(in), optional :: check_shift_consistency
            !! Check consitency of shift array
        
        logical :: do_check_uniqueness, do_check_shift_consistency
        
        ! Set defaults for optional parameters
        do_check_uniqueness = .true.
        if (present(check_uniqueness)) do_check_uniqueness = check_uniqueness
        
        do_check_shift_consistency = .true.
        if (present(check_shift_consistency)) do_check_shift_consistency = check_shift_consistency
        
        call set_ok(ierr)
                
        ! 1. Check basic structure
        call validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                                   gene_to_fam, expression_vectors, family_centroids, &
                                   shift_vectors, ierr)
        if (.not. is_ok(ierr)) return
        
        ! 2. Check gene-to-family mapping
        call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        if (.not. is_ok(ierr)) return
        
        ! 3. Check expression data
        call validate_expression_data(expression_vectors, .true., ierr)  ! Check for non-negative
        if (.not. is_ok(ierr)) return
        
        ! 4. Check family centroids
        if (n_families > 0) then
            call validate_family_centroids(family_centroids, ierr)
            if (.not. is_ok(ierr)) return
        else 
            if(DEBUG) write(*,*) 'Warning: No families defined, skipping family centroid checks.'
        end if
        
        ! 5. Check shift vectors
        if (do_check_shift_consistency) then 
            call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                  gene_to_fam, n_samples, ierr)
            if (.not. is_ok(ierr)) return
        end if

        ! 6. Check uniqueness (optional, can be slow for large datasets)
        if (do_check_uniqueness) then
            call validate_string_array_uniqueness(gene_ids, ierr)
            if (.not. is_ok(ierr)) return
            
            if (n_families > 0) then
                call validate_string_array_uniqueness(gene_family_ids, ierr)
                if (.not. is_ok(ierr)) return
            end if
        end if
        
        if(DEBUG) write(*,*) 'All data validation checks passed!'
        
    end subroutine validate_all_data
    
end module tox_data_validation

!> R binding for gene to family mapping validating
subroutine validate_gene_to_family_mapping_R(gene_to_fam, n_genes, n_families, ierr)
    use tox_data_validation, only: validate_gene_to_family_mapping
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: n_families
        !! Number of families
    integer(int32), intent(in) :: gene_to_fam(n_genes)
        !! Gene to family mapping

    integer(int32), intent(out) :: ierr
        !! Error code
    call set_ok(ierr)
    call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
end subroutine validate_gene_to_family_mapping_R

!> R Binding to validate expression data
subroutine validate_expression_data_R(expression_vectors, n_genes, n_samples, check_non_negative, ierr)
    use tox_data_validation, only: validate_expression_data
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    logical, intent(in) :: check_non_negative
        !! Defines if expression data should be checked for negative values
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    real(real64), intent(in) :: expression_vectors(n_genes, n_samples)
        !! Expression vectors array
    integer(int32), intent(out) :: ierr
        !! Error code
    call set_ok(ierr)
    call validate_expression_data(expression_vectors, check_non_negative, ierr)
end subroutine validate_expression_data_R

!> R Binding to validate family centroids
subroutine validate_family_centroids_R(family_centroids, n_families, n_samples, ierr)
    use tox_data_validation, only : validate_family_centroids
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_families
        !! Number of families
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    real(real64), intent(in) :: family_centroids(n_samples, n_families)
        !! Family centroids array
    integer(int32), intent(out) :: ierr
        !! Error code
    call set_ok(ierr)
    call validate_family_centroids(family_centroids, ierr)
end subroutine validate_family_centroids_R

!> R Binding to validate shift vectors
subroutine validate_shift_vectors_R(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_genes, n_samples, n_families, ierr)
    use tox_data_validation, only: validate_shift_vectors
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: n_families
        !! Number of families
    real(real64), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! shift vectors array
    real(real64), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Expression vectors array
    real(real64), intent(in) :: family_centroids(n_samples, n_families)
        !! family centroids array
    integer(int32), intent(in) :: gene_to_fam(n_genes)
        !! gene to family mapping
    integer(int32), intent(out) :: ierr
        !! Error code
    call set_ok(ierr)
    
    call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples, ierr)
end subroutine validate_shift_vectors_R

!> R Binding to validate string array uniqueness
subroutine validate_string_array_uniqueness_R(str_arr, str_len, n_strings, ierr)
    use tox_data_validation, only: validate_string_array_uniqueness
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: str_len
        !! Length of the strings
    integer(int32), intent(in) :: n_strings
        !! Number of strings
    character(kind=c_char, len=1), intent(in) :: str_arr(str_len, n_strings)
        !! strings array
    integer(int32), intent(out) :: ierr
        !! Error code

    character(len=:), allocatable :: temp_str_arr(:)
    character(len=:), allocatable :: temp_str
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=str_len) :: temp_str_arr(n_strings), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if

    do i = 1, n_strings
        call c_char_1d_as_string(str_arr(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        temp_str_arr(i) = temp_str
    end do

    call validate_string_array_uniqueness(temp_str_arr, ierr)
end subroutine validate_string_array_uniqueness_R

!> R Binding to validate data structure
subroutine validate_data_structure_R(n_genes, n_families, n_samples, &
                                     gene_ids_raw, gene_ids_len, &
                                     gene_family_ids_raw, fam_len, &
                                     gene_to_fam, expression_vectors, family_centroids, &
                                     shift_vectors, ierr)
    use tox_data_validation, only: validate_data_structure
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: n_families
        !! Number of families
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    integer(int32), intent(in) :: gene_ids_len
        !! gene ids length    
    character(kind=c_char, len=1), intent(in) :: gene_ids_raw(gene_ids_len, n_genes)
        !! gene ids array
    integer(int32), intent(in) :: fam_len
        !! Length of the family ids
    character(kind=c_char, len=1), intent(in) :: gene_family_ids_raw(fam_len, n_genes)
        !! family ids array
    integer(int32), intent(in) :: gene_to_fam(n_genes)
        !! gene to family array
    real(real64), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Expression vectors array
    real(real64), intent(in) :: family_centroids(n_samples, n_families)
        !! family centroids array
    real(real64), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! shift vectors array
    integer(int32), intent(out) :: ierr
        !! Error code
    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=gene_ids_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_genes), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if
    do i = 1, n_genes
        call c_char_1d_as_string(gene_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_ids(i) = temp_str
        
        call c_char_1d_as_string(gene_family_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_family_ids(i) = temp_str
    end do

    call validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                                 gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_data_structure_R

!> R Binding to validate all data
subroutine validate_all_data_R(n_genes, n_families, n_samples, &
                               gene_ids_raw, gene_len, &
                               gene_family_ids_raw, fam_len, &
                               gene_to_fam, expression_vectors, family_centroids, &
                               shift_vectors, ierr)
    use tox_data_validation, only: validate_all_data
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: n_families
        !! Number of families
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    integer(int32), intent(in) :: gene_len
        !! length of the gene ids    
    character(kind=c_char, len=1), intent(in) :: gene_ids_raw(gene_len, n_genes)
        !! Gene ids array
    integer(int32), intent(in) :: fam_len
        !! Length of the family indices
    character(kind=c_char, len=1), intent(in) :: gene_family_ids_raw(fam_len, n_families)
        !! Family ids array
    integer(int32), intent(in) :: gene_to_fam(n_genes)
        !! genes to family mapping
    real(real64), intent(in) :: expression_vectors(n_samples, n_genes)
        !! expression vectors
    real(real64), intent(in) :: family_centroids(n_samples, n_families)
        !! family centroids
    real(real64), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! Shift vectors
    integer(int32), intent(out) :: ierr
        !! error code
    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=gene_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_families), stat=ios)
    do i = 1, n_genes
        call c_char_1d_as_string(gene_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_ids(i) = temp_str
    end do

    do i=1, n_families
        call c_char_1d_as_string(gene_family_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_family_ids(i) = temp_str
    end do

    call validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                           gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_all_data_R

! ---- C bindings for validation routines ----

!> C binding to validate gene to family mapping
subroutine validate_gene_to_family_mapping_C(gene_to_fam, n_genes, n_families, ierr) bind(C, name="validate_gene_to_family_mapping_C")
    use iso_c_binding, only: c_int
    use tox_data_validation, only: validate_gene_to_family_mapping
    use tox_errors, only: set_ok
    implicit none
    integer(c_int), intent(in), value :: n_genes
        !! Number of genes 
    integer(c_int), intent(in) :: gene_to_fam(n_genes)
        !! Pointer to gene to family array
    integer(c_int), intent(in), value :: n_families
        !! Numbero of families
    integer(c_int), intent(out) :: ierr
        !! Error code
        
    call set_ok(ierr)
    call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
end subroutine validate_gene_to_family_mapping_C

!> C Binding to validate expression data
subroutine validate_expression_data_C(expression_vectors, n_genes, n_samples, check_non_negative, ierr) bind(C, name="validate_expression_data_C")
    use iso_c_binding, only: c_int, c_double
    use tox_data_validation, only: validate_expression_data
    use tox_errors, only: set_ok
    use tox_conversions, only: c_int_as_logical
    implicit none
    integer(c_int), intent(in), value :: n_genes
        !! Number of genes
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples
    real(c_double), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors array
    integer(c_int), intent(in), value :: check_non_negative
        !! Define if expression data should be checked for negative values
    integer(c_int), intent(out) :: ierr
        !! Error code

    logical :: l_check
    call set_ok(ierr)
    call c_int_as_logical(check_non_negative, l_check)
    call validate_expression_data(expression_vectors, l_check, ierr)
end subroutine validate_expression_data_C

!> C binding to validate family centroids
subroutine validate_family_centroids_C(family_centroids, n_families, n_samples, ierr) bind(C, name="validate_family_centroids_C")
    use iso_c_binding, only: c_int, c_double
    use tox_data_validation, only: validate_family_centroids
    use tox_errors, only: set_ok
    implicit none
    integer(c_int), intent(in), value :: n_families
        !! Number of families
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples
    real(c_double), intent(in) :: family_centroids(n_samples, n_families)
        !! Pointer to family_centroids array
    integer(c_int), intent(out) :: ierr
        !! Error code

    call set_ok(ierr)
    call validate_family_centroids(family_centroids, ierr)
end subroutine validate_family_centroids_C

!> C Binding to validate shift vectors
subroutine validate_shift_vectors_C(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_genes, n_samples, n_families, ierr) bind(C, name="validate_shift_vectors_C")
    use iso_c_binding, only: c_int, c_double
    use tox_data_validation, only: validate_shift_vectors
    use tox_errors, only: set_ok
    implicit none
    integer(c_int), intent(in), value :: n_genes
        !! Number of genes
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples
    integer(c_int), intent(in), value :: n_families
        !! Number of families
    real(c_double), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! Pointer to shift vectors array
    real(c_double), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors array
    real(c_double), intent(in) :: family_centroids(n_samples, n_families)
        !! Pointer to family centroids array
    integer(c_int), intent(in) :: gene_to_fam(n_genes)
        !! Pointer to gene to family array
    integer(c_int), intent(out) :: ierr
        !! Error code

    call set_ok(ierr)

    call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples, ierr)
end subroutine validate_shift_vectors_C

!> C Binding to validate string array uniqueness
subroutine validate_string_array_uniqueness_C(str_arr, str_len, n_strings, ierr) bind(C, name="validate_string_array_uniqueness_C")
    use iso_c_binding, only: c_int, c_char
    use tox_data_validation, only: validate_string_array_uniqueness
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    implicit none
    integer(c_int), intent(in), value :: str_len
        !! Length of the strings
    integer(c_int), intent(in), value :: n_strings
        !! Number of strings
    character(kind=c_char, len=1), intent(in) :: str_arr(str_len, n_strings)
        !! Pointer to string array
    integer(c_int), intent(out) :: ierr 
        !! Error code

    character(len=:), allocatable :: temp_str_arr(:)
    character(len=:), allocatable :: temp_str
    integer :: i, ios
    
    call set_ok(ierr)
    allocate(character(len=str_len) :: temp_str_arr(n_strings), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if
    do i = 1, n_strings
        call c_char_1d_as_string(str_arr(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        temp_str_arr(i) = temp_str
    end do
    call validate_string_array_uniqueness(temp_str_arr, ierr)
end subroutine validate_string_array_uniqueness_C

!> C Binding to validate data structure
subroutine validate_data_structure_C(n_genes, n_families, n_samples, &
                                     gene_ids_raw, gene_ids_len, &
                                     gene_family_ids_raw, fam_len, &
                                     gene_to_fam, expression_vectors, family_centroids, &
                                     shift_vectors, ierr) bind(C, name="validate_data_structure_C")
    use iso_c_binding, only: c_int, c_double, c_char
    use tox_data_validation, only: validate_data_structure
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    implicit none
    integer(c_int), intent(in), value :: n_genes
        !! Number of genes
    integer(c_int), intent(in), value :: n_families
        !! Number of families
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples
    integer(c_int), intent(in), value :: gene_ids_len   
        !! Length of gene ids
    integer(c_int), intent(in), value :: fam_len
        !! Length of the family indices        
    character(kind=c_char, len=1), intent(in) :: gene_ids_raw(gene_ids_len, n_genes)
        !! Pointer to gene ids array
    character(kind=c_char, len=1), intent(in) :: gene_family_ids_raw(fam_len, n_families)
        !! Pointer to family ids array
    integer(c_int), intent(in) :: gene_to_fam(n_genes)
        !! Pointer to gene to family array
    real(c_double), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors
    real(c_double), intent(in) :: family_centroids(n_samples, n_families)
        !! Pointer to family centroids
    real(c_double), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! Pointer to shift vectors
    integer(c_int), intent(out) :: ierr
        !! Error code

    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer :: i, ios
    
    call set_ok(ierr)

    allocate(character(len=gene_ids_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_families), stat=ios)

    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if

    do i = 1, n_genes
        call c_char_1d_as_string(gene_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_ids(i) = temp_str
    end do 
    do i = 1, n_families
        call c_char_1d_as_string(gene_family_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_family_ids(i) = temp_str
    end do

    call validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                                 gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_data_structure_C

!> C Binding to validate all data
subroutine validate_all_data_C(n_genes, n_families, n_samples, &
                               gene_ids_raw, gene_len, &
                               gene_family_ids_raw, fam_len, &
                               gene_to_fam, expression_vectors, family_centroids, &
                               shift_vectors, ierr) bind(C, name="validate_all_data_C")
    use iso_c_binding, only: c_int, c_double, c_char
    use tox_data_validation, only: validate_all_data
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use tox_conversions, only: c_char_1d_as_string
    implicit none
    integer(c_int), intent(in), value :: n_genes
        !! Number of genes
    integer(c_int), intent(in), value :: n_families
        !! Number of families
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples
    integer(c_int), intent(in), value :: gene_len
        !! Length of the gene ids    
    character(kind=c_char, len=1), intent(in) :: gene_ids_raw(gene_len, n_genes)
        !! Pointer to gene ids array
    integer(c_int), intent(in), value :: fam_len
        !! Length of the family ids    
    character(kind=c_char, len=1), intent(in) :: gene_family_ids_raw(fam_len, n_families)
        !! pointer to family ids array
    integer(c_int), intent(in) :: gene_to_fam(n_genes)
        !! Pointer to gene to family mapping
    real(c_double), intent(in) :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors
    real(c_double), intent(in) :: family_centroids(n_samples, n_families)
        !! Pointer to family centroids
    real(c_double), intent(in) :: shift_vectors(2*n_samples, n_genes)
        !! Pointer to shift vectors
    integer(c_int), intent(out) :: ierr
        !! Error code

    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer :: i, ios
    
    call set_ok(ierr)

    allocate(character(len=gene_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_families), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if
    do i = 1, n_genes
        call c_char_1d_as_string(gene_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_ids(i) = temp_str
    end do
    do i=1, n_families
        call c_char_1d_as_string(gene_family_ids_raw(:,i), temp_str, ierr)
        if (.not. is_ok(ierr)) return
        gene_family_ids(i) = temp_str
    end do
    call validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                           gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_all_data_C