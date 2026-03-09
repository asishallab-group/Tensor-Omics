! filepath: test/mod_test_tox_data.f90
!> Unit test suite for expression readers and data processing routines.
module mod_test_tox_data
  use asserts
  use iso_fortran_env, only: real64, int32
  use tox_data_tools
  use tox_data_validation
  use tox_data_accessors
  use tox_data_read_write
  use f42_array_utils, only: get_array_metadata
  use tox_gene_centroids
  use tox_shift_vectors
  use tox_errors
  use tox_archive, only: save_tox_data, read_tox_data, create_zip_archive, extract_zip_archive, delete_file
  use f42_deserialize_int
  use f42_deserialize_real
  use f42_deserialize_char
  use f42_serialize_char
  use f42_serialize_int
  use f42_serialize_real
  use mod_test_suite, only: test_case
  implicit none
  public


  ! Global test data
  integer(int32), allocatable :: gene_to_fam(:)
  character(len=256), allocatable :: gene_family_ids(:), gene_ids(:)
  real(real64), allocatable :: kallisto_expr(:,:), kallisto_expr_verify(:,:)
  integer(int32) :: n_genes, n_families, total_samples
  real(real64), allocatable :: family_centroids(:,:)
  real(real64), allocatable :: shift_vectors(:,:)  ! Added for shift vectors

contains

  !> Get array of all available tests.
  function get_all_tests_tox_data() result(all_tests)
    type(test_case),allocatable:: all_tests(:)
    allocate(all_tests(23))
    all_tests(1) = test_case("test_read_gene_ids", test_read_gene_ids)
    all_tests(2) = test_case("test_read_expression_data", test_read_expression_data)
    all_tests(3) = test_case("test_read_family_mapping", test_read_family_mapping)
    all_tests(4) = test_case("test_validate_data", test_validate_data)
    all_tests(5) = test_case("test_compute_centroids", test_compute_centroids)
    all_tests(6) = test_case("test_write_read_expression_data", test_write_read_expression_data)
    all_tests(7) = test_case("test_compute_shift_vectors", test_compute_shift_vectors)
    all_tests(8) = test_case("test_write_read_shift_vectors", test_write_read_shift_vectors)
    all_tests(9) = test_case("test_read_write_gene_ids", test_read_write_gene_ids)
    all_tests(10) = test_case("test_read_write_gene_to_fam", test_read_write_gene_to_fam)
    all_tests(11) = test_case("test_read_write_family_ids", test_read_write_family_ids)
    all_tests(12) = test_case("test_read_write_family_centroids", test_read_write_centroids)
    all_tests(13) = test_case("test_data_accessors", test_data_accessors)
    all_tests(14) = test_case("test_archive", test_archive)
    all_tests(15) = test_case("test_hashing", test_hashing)
    all_tests(16) = test_case("test_validate_duplicate_gene_ids", test_validate_duplicate_gene_ids)
    all_tests(17) = test_case("test_validate_nan_direct", test_validate_nan_direct)
    all_tests(18) = test_case("test_validate_inf_direct", test_validate_inf_direct)
    all_tests(19) = test_case("test_validate_negative_expression", test_validate_negative_expression)
    all_tests(20) = test_case("test_validate_invalid_gene_to_fam", test_validate_invalid_gene_to_fam)
    all_tests(21) = test_case("test_validate_empty_strings", test_validate_empty_strings)
    all_tests(22) = test_case("test_validate_dimension_mismatch", test_validate_dimension_mismatch)
    all_tests(23) = test_case("test_manual_archive", test_manual_archive)
  end function get_all_tests_tox_data

   
  

  !> Setup global test data
  subroutine setup_global_data()
    character(len=256), allocatable :: expr_file(:)
    integer(int32) :: ierr, i, n_genes_kept
    logical, allocatable :: ortholog_mask(:), unassigned_mask(:)
    integer(int32), allocatable :: selected_indices(:), value_cols(:)

    ! Initialize file lists
    allocate(expr_file(1))
    allocate(value_cols(67))

    do i = 2, 68
      value_cols(i-1) = i
    end do
    
    expr_file = ['material/kallisto_sex_data_no_na.tsv']
    ! Calculate total samples
    total_samples = 67
    n_genes = 88327  ! Original number of genes
    n_families = 15512

    ! Allocate arrays with original size
    allocate(gene_ids(n_genes))
    allocate(kallisto_expr(total_samples, n_genes))
    allocate(gene_family_ids(n_families))
    allocate(gene_to_fam(n_genes))
    allocate(unassigned_mask(n_genes))

    ! Read gene IDs
    call read_gene_ids_from_tsv_file(expr_file(1), gene_ids, 1, 1, ierr)
    call assert_equal_int(ierr, 0, "Reading gene IDs should succeed")

    ! write(*,*) 'First 10 gene IDs:'
    ! do i = 1, min(10, n_genes)
    !   write(*,*) trim(gene_ids(i))
    ! end do

    ! Read expression data
    kallisto_expr = 0.0_real64
    call read_expression_vectors_tsv(expr_file, gene_ids, kallisto_expr, &
                          1, 1, value_cols, 1, ierr, char(9))

    ! Read family mapping
    call read_orthofinder_file('material/Orthogroups.tsv', gene_ids, gene_family_ids, gene_to_fam, ierr)
    call assert_equal_int(ierr, ERR_OK, "Reading family file should succeed")
    
    ! Print first 10 family IDs
    ! write(*,*) 'First 10 family IDs:'
    ! do i = 1, min(10, n_families)
    !   write(*,*) trim(gene_family_ids(i))
    ! end do
    ! Filter out genes without family assignments
    call get_unassigned_mask(gene_to_fam, unassigned_mask, n_genes_kept)
    call apply_unassigned_mask(gene_ids, kallisto_expr, gene_to_fam, unassigned_mask, n_genes_kept, ierr)

    call assert_equal_int(ierr, ERR_OK, "Filtering unassigned genes should succeed")
    n_genes = n_genes_kept  ! Update n_genes to reflect the filtered count

    ! Compute centroids
    allocate(family_centroids(total_samples, n_families))
    
    ! Create explicit arrays instead of array constructors
    allocate(ortholog_mask(n_genes))
    allocate(selected_indices(n_genes))
    
    do i = 1, n_genes
      ortholog_mask(i) = .true.
      selected_indices(i) = i
    end do

    ! write(*,*) 'Size of family_centroids: ', size(family_centroids, 1), size(family_centroids, 2)
    
    call group_centroid(kallisto_expr, total_samples, n_genes, gene_to_fam, &
                      n_families, family_centroids, GROUP_ALL, selected_indices, ierr, ortholog_mask)
    call assert_equal_int(ierr, ERR_OK, "Computing centroids should succeed")
    
    ! Compute shift vectors
    allocate(shift_vectors(2*total_samples, n_genes))
    call compute_shift_vector_field(total_samples, n_genes, n_families, kallisto_expr, &
                                  family_centroids, gene_to_fam, shift_vectors, ierr)
    call assert_equal_int(ierr, ERR_OK, "Computing shift vectors should succeed")
    
    ! Clean up temporary arrays
    deallocate(ortholog_mask, selected_indices)
  end subroutine setup_global_data

 

  !> Test reading gene IDs
  subroutine test_read_gene_ids()
    if (.not. allocated(gene_ids)) then
      call setup_global_data()
    end if
    call assert_true(allocated(gene_ids), "Gene IDs should be allocated")
    call assert_equal_int(size(gene_ids), n_genes, "Number of gene IDs should match")
    call assert_true(len_trim(gene_ids(1)) > 0, "First gene ID should not be empty")
  end subroutine test_read_gene_ids

  !> Test reading expression data
  subroutine test_read_expression_data()
    character(len=256), allocatable :: gene_ids_false_inputs(:)
    real(real64), allocatable :: expr_vecs_false_inputs(:,:)
    character(len=64), allocatable :: inf_file(:), nan_file(:), missing_col_file(:), invalid_struct_file(:), no_genes_file(:), empty_file(:), mixed_seperators(:)
    character(len=64), allocatable :: csv_file(:)
    integer(int32) :: ierr

    call set_ok(ierr)

    allocate(inf_file(1))
    allocate(nan_file(1))
    allocate(missing_col_file(1))
    allocate(gene_ids_false_inputs(10))
    allocate(expr_vecs_false_inputs(6, n_genes))
    allocate(invalid_struct_file(1))
    allocate(no_genes_file(1))
    allocate(empty_file(1))
    allocate(mixed_seperators(1))
    allocate(csv_file(1))

    call assert_true(allocated(kallisto_expr), "Expression data should be allocated")
    call assert_equal_int(size(kallisto_expr, 1), total_samples, "Number of samples should match")
    call assert_equal_int(size(kallisto_expr, 2), n_genes, "Number of genes should match")
    call assert_true(all(kallisto_expr >= 0.0_real64), "All expression values should be non-negative")

    call read_gene_ids_from_tsv_file('test/test_files/kallisto_dup_gene_ids.tsv', gene_ids_false_inputs, 1, 1, ierr)
    call assert_equal_int(ierr, ERR_OK, "Error while reading duplicate gene ids")

    call validate_string_array_uniqueness(gene_ids_false_inputs, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Duplicate gene IDs should be detected")

    inf_file = [ &
      'test/test_files/kallisto_Inf.tsv']

    nan_file = [ &
      'test/test_files/kallisto_NaN.tsv']
    
    invalid_struct_file = [ &
      'test/test_files/kallisto_invalid_strucutre.tsv']

    no_genes_file = [ &
      'test/test_files/kallisto_no_gene_ids.tsv']
    
    empty_file = [ &
      'test/test_files/kallisto_empty_file.tsv']

    mixed_seperators = [ &
      'test/test_files/kallisto_mixed_seperators.tsv']
    
    csv_file = [ &
      'test/test_files/kallisto_csv_values.csv']

    call read_expression_vectors_tsv(inf_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [2,3,4,5,6,7], 1, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Error while reading expression vectors, should get invalid input for Inf in expression data")
    call set_ok(ierr)

    call read_expression_vectors_tsv(nan_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [2,3,4,5,6,7], 1, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Error while reading expression vectors, should get invalid input for NaN in expression data")
    call set_ok(ierr)

    call read_expression_vectors_tsv(invalid_struct_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [2,3,4,5,6,7], 1, ierr)
    call assert_equal_int(ierr, ERR_READ_DATA, "Should throw error for invalid structure")
    call set_ok(ierr)

    call read_expression_vectors_tsv(no_genes_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [1,2,3,4,5,6], 1, ierr)
    call assert_equal_int(ierr, ERR_OK, "Should print warnings for missing genes")
    call set_ok(ierr)

    call read_expression_vectors_tsv(empty_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [1,2,3,4,5,6], 1, ierr)
    call assert_equal_int(ierr, ERR_READ_DATA, "should throw error for empty file")
    call set_ok(ierr)

    call read_expression_vectors_tsv(mixed_seperators, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [2,3,4,5,6,7], 1, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Should throw error for file with mixed seperators")
    call set_ok(ierr)

    call read_expression_vectors_tsv(csv_file, gene_ids_false_inputs, expr_vecs_false_inputs, 1, 1, [2,3,4,5,6,7], 1, ierr, ',')
    call assert_equal_int(ierr, ERR_OK, "Should read csv file without error")

  end subroutine test_read_expression_data

  !> Test reading family mapping
  subroutine test_read_family_mapping()
    call assert_true(allocated(gene_to_fam), "Gene to family mapping should be allocated")
    call assert_true(allocated(gene_family_ids), "Family IDs should be allocated")
    call assert_equal_int(size(gene_to_fam), n_genes, "Gene to family mapping size should match")
    call assert_equal_int(size(gene_family_ids), n_families, "Number of family IDs should match")
  end subroutine test_read_family_mapping

  !> Test data validation
  subroutine test_validate_data()
    integer(int32) :: ierr

    call validate_data_structure(n_genes, n_families, total_samples, gene_ids, gene_family_ids, gene_to_fam, &
                                kallisto_expr, family_centroids, shift_vectors, ierr)
    call assert_equal_int(ierr, ERR_OK, "Data structure could not be validated")

    call validate_empty_strings(gene_ids, "Gene IDs", ierr)
    call assert_equal_int(ierr, ERR_OK, "Gene IDs could not be validated - contains empty strings")

    call validate_empty_strings(gene_family_ids, "Family IDs", ierr)
    call assert_equal_int(ierr, ERR_OK, "Family Ids contain empty entries")

    call validate_expression_data(kallisto_expr, .true., ierr)
    call assert_equal_int(ierr, ERR_OK, "Expression data could not be validated")

    call validate_all_data(n_genes, n_families, total_samples, gene_ids, gene_family_ids, &
                           gene_to_fam, kallisto_expr, family_centroids, shift_vectors, ierr, .true., .true.)
    call assert_equal_int(ierr, ERR_OK, "Data could not be validated")
  end subroutine test_validate_data

  !> Test centroid computation
  subroutine test_compute_centroids()
    integer(int32) :: i, j, n_family_genes, current_idx
    integer(int32), allocatable :: indice_set(:)
    real(real64), allocatable :: test_centroid(:)
    real(real64) :: sample_sum
    
    call assert_true(allocated(family_centroids), "Family centroids should be allocated")
    call assert_equal_int(size(family_centroids, 1), total_samples, "Centroid samples should match")
    call assert_equal_int(size(family_centroids, 2), n_families, "Number of centroid families should match")
    
    ! Test specific family (family 2)
    n_family_genes = count(gene_to_fam == 2)
    allocate(indice_set(n_family_genes))
    allocate(test_centroid(total_samples))
    
    current_idx = 1
    do i = 1, n_genes
      if (gene_to_fam(i) == 2) then
        indice_set(current_idx) = i
        current_idx = current_idx + 1
      end if
    end do

    do j = 1, total_samples
      sample_sum = 0.0_real64
      do i = 1, n_family_genes
        sample_sum = sample_sum + kallisto_expr(j, indice_set(i))
      end do
      test_centroid(j) = sample_sum / n_family_genes
    end do

    ! Fixed call to match the assert_equal_array_real interface
    call assert_equal_array_real(test_centroid, family_centroids(:, 2), &
                                total_samples, 1e-12_real64, &
                                "Family 2 centroid should match manual calculation")
  end subroutine test_compute_centroids

  !> Test binary write/read operations
  subroutine test_write_read_expression_data()
    integer(int32) :: ierr, ndims, dims(2)
    integer(int32) :: total_elements
    
    call save_expression_vectors(kallisto_expr, 'test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving expression data should succeed")
    
    call get_array_metadata('test_kallisto_data.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Reading metadata should succeed")
    call assert_equal_int(ndims, 2, "Array should have 2 dimensions")
    call assert_equal_int(dims(1), total_samples, "First dimension should match sample count")
    call assert_equal_int(dims(2), n_genes, "Second dimension should match gene count")
    
    allocate(kallisto_expr_verify(dims(1), dims(2)))
    call load_expression_vectors(kallisto_expr_verify, 'test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading expression data should succeed")
    
    ! Fixed call to match the assert_equal_array_real interface
    total_elements = size(kallisto_expr)
    call assert_equal_array_real(kallisto_expr, &
                                kallisto_expr_verify, &
                                total_elements, 1e-12_real64, &
                                "Original and loaded data should be identical")
    
    ! Clean up
    deallocate(kallisto_expr_verify)
  end subroutine test_write_read_expression_data

  !> Test computation of shift vectors
  subroutine test_compute_shift_vectors()
    integer(int32) :: i
    
    call assert_true(allocated(shift_vectors), "Shift vectors should be allocated")
    call assert_equal_int(size(shift_vectors, 1), 2*total_samples, "Shift vectors should have 2*d rows")
    call assert_equal_int(size(shift_vectors, 2), n_genes, "Shift vectors should have n_genes columns")
    
    ! Verify centroids are correctly placed in first d rows
    do i = 1, n_genes
      if (gene_to_fam(i) >= 1 .and. gene_to_fam(i) <= n_families) then
        call assert_equal_array_real(shift_vectors(1:total_samples, i), &
                                    family_centroids(:, gene_to_fam(i)), &
                                    total_samples, 1e-12_real64, &
                                    "Centroid should match family centroid")
      end if
    end do
    
    ! Verify shift vectors are correctly computed
    do i = 1, n_genes
      if (gene_to_fam(i) >= 1 .and. gene_to_fam(i) <= n_families) then
        call assert_equal_array_real(shift_vectors(total_samples+1:2*total_samples, i), &
                                    kallisto_expr(:, i) - family_centroids(:, gene_to_fam(i)), &
                                    total_samples, 1e-12_real64, &
                                    "Shift vector should be expression minus centroid")
      end if
    end do
  end subroutine test_compute_shift_vectors

  !> Test writing and reading shift vectors
  subroutine test_write_read_shift_vectors()
    integer(int32) :: ierr, ndims, dims(2)
    integer(int32) :: total_elements
    real(real64), allocatable :: shift_vectors_loaded(:,:)
    
    call save_expression_vectors(shift_vectors, 'test_shift_vectors.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving shift vectors should succeed")
    
    call get_array_metadata('test_shift_vectors.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Reading metadata should succeed")
    call assert_equal_int(ndims, 2, "Array should have 2 dimensions")
    call assert_equal_int(dims(1), 2*total_samples, "First dimension should match 2*d")
    call assert_equal_int(dims(2), n_genes, "Second dimension should match gene count")
    
    allocate(shift_vectors_loaded(dims(1), dims(2)))
    call load_expression_vectors(shift_vectors_loaded, 'test_shift_vectors.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading shift vectors should succeed")
    
    ! Verify loaded data matches original
    total_elements = size(shift_vectors)
    call assert_equal_array_real(shift_vectors, &
                                shift_vectors_loaded, &
                                total_elements, 1e-12_real64, &
                                "Original and loaded shift vectors should be identical")
    
    ! Clean up
    deallocate(shift_vectors_loaded)
  end subroutine test_write_read_shift_vectors

  subroutine test_read_write_gene_ids()
    integer(int32) :: ierr, ndims, dims(1)
    character(len=256), allocatable :: loaded_gene_ids(:)
    call save_gene_ids(gene_ids, 'test_gene_ids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving gene IDs should succeed")

    call get_array_metadata('test_gene_ids.bin', dims, 1, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Getting metadata should succeed")

    call assert_equal_int(ndims, 1, "Gene IDs should be 1D array")
    allocate(loaded_gene_ids(dims(1)))
    call load_gene_ids(loaded_gene_ids, 'test_gene_ids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading gene IDs should succeed")

    call assert_equal_array_char(loaded_gene_ids, gene_ids, 128, dims(1), &
                                "Loaded gene IDs should match original")
  end subroutine test_read_write_gene_ids

  subroutine test_read_write_gene_to_fam()
    integer(int32) :: ierr, ndims, dims(1)
    integer(int32), allocatable :: loaded_gene_to_fam(:)
    call save_gene_to_family(gene_to_fam, 'test_gene_to_fam.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving gene to family mapping should succeed")
    call get_array_metadata('test_gene_to_fam.bin', dims, 1, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Getting metadata should succeed")
    call assert_equal_int(ndims, 1, "Gene to family mapping should be 1D array")
    allocate(loaded_gene_to_fam(dims(1)))
    call load_gene_to_family(loaded_gene_to_fam, 'test_gene_to_fam.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading gene to family mapping should succeed")
    call assert_equal_array_int(loaded_gene_to_fam, gene_to_fam, dims(1), &
                               "Loaded gene to family mapping should match original")
  end subroutine test_read_write_gene_to_fam

  subroutine test_read_write_family_ids()
    integer(int32) :: ierr, ndims, dims(1)
    character(len=256), allocatable :: loaded_family_ids(:)
    call save_family_ids(gene_family_ids, 'test_family_ids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving family IDs should succeed")
    call get_array_metadata('test_family_ids.bin', dims, 1, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Getting metadata should succeed")
    call assert_equal_int(ndims, 1, "Family IDs should be 1D array")
    allocate(loaded_family_ids(dims(1)))
    call load_family_ids(loaded_family_ids, 'test_family_ids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading family IDs should succeed")
    call assert_equal_array_char(loaded_family_ids, gene_family_ids, 128, dims(1), &
                                "Loaded family IDs should match original")
  end subroutine test_read_write_family_ids

  subroutine test_read_write_centroids()
    integer(int32) :: ierr, ndims, dims(2)
    real(real64), allocatable :: loaded_centroids(:,:)
    call save_family_centroids(family_centroids, 'test_family_centroids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Saving family centroids should succeed")
    call get_array_metadata('test_family_centroids.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Getting metadata should succeed")
    call assert_equal_int(ndims, 2, "Family centroids should be 2D array")
    allocate(loaded_centroids(dims(1), dims(2)))
    call load_family_centroids(loaded_centroids, 'test_family_centroids.bin', ierr)
    call assert_equal_int(ierr, ERR_OK, "Loading family centroids should succeed")
    call assert_equal_array_real(loaded_centroids, &
                                family_centroids, &
                                size(family_centroids), 1e-12_real64, &
                                "Loaded family centroids should match original")
  end subroutine test_read_write_centroids

  subroutine test_data_accessors()
    real(real64), allocatable :: returned_centroid(:)
    integer(int32) :: family_idx, ierr, gene_idx

    allocate(returned_centroid(total_samples))

    gene_idx = get_gene_index(gene_ids, 'NP_001000001.1')
    call assert_equal_int(gene_idx, 2, "Gene index for NP_001000001.1 should be 2")
    call get_family_for_gene_index(gene_idx, gene_to_fam, family_idx, ierr)
    call assert_equal_int(ierr, ERR_OK, "Getting family index should succeed")
    call assert_equal_int(family_idx, 48, "Family index for gene index 2")
    call assert_string_equal(gene_family_ids(family_idx), 'OG0000047', &
                        "Family ID for family index 48 should be OG0000047")

    family_idx = get_family_index(gene_family_ids, 'OG0000001')
    call assert_equal_int(family_idx, 2, "Family index for OG0000001 should be 2")
    call get_family_centroid(family_idx, family_centroids, returned_centroid, ierr)

    call assert_equal_int(ierr, ERR_OK, "Getting family centroid should succeed")
    call assert_equal_array_real(returned_centroid, family_centroids(:, family_idx), &
                               size(family_centroids, 1), 1e-12_real64, &
                               "Retrieved centroid should match original")
  end subroutine test_data_accessors

  subroutine test_archive()
    use iso_fortran_env, only: real64, int32
    implicit none
    
    ! Declare verification arrays
    character(len=:), allocatable :: gene_ids_verify(:)
    real(real64), allocatable :: kallisto_verify(:,:)
    integer(int32), allocatable :: gene_to_fam_verify(:)
    character(len=:), allocatable :: gene_family_ids_verify(:)
    real(real64), allocatable :: family_centroids_verify(:,:)
    real(real64), allocatable :: shift_vectors_verify(:,:)
    
    integer(int32) :: ierr
    
    ! Test 1: Save and read all data
    ! print *, "Test 1: Saving and reading all data"
    call save_tox_data("test_archive_1_f.zip", ierr, &
                      gene_ids=gene_ids, gene_ids_file="gene_ids_v1.bin", &
                      expression=kallisto_expr, expression_file="kallisto_v1.bin", &
                      gene_to_family=gene_to_fam, gene_to_family_file="gene_to_fam.bin", &
                      family_ids=gene_family_ids, family_ids_file="family_ids.bin", &
                      family_centroids=family_centroids, family_centroids_file="family_centroids.bin", &
                      shift_vectors=shift_vectors, shift_vectors_file="shift_vectors.bin")
    
    call assert_equal_int(ierr, ERR_OK, "Error saving archive")
    
    call read_tox_data("test_archive_1_f.zip", ierr, &
                      gene_ids=gene_ids_verify, &
                      expression=kallisto_verify, &
                      gene_to_family=gene_to_fam_verify, &
                      family_ids=gene_family_ids_verify, &
                      family_centroids=family_centroids_verify, &
                      shift_vectors=shift_vectors_verify)
    
    call assert_equal_int(ierr, ERR_OK, "Error reading archive")

    call assert_true(allocated(gene_ids_verify), "gene_ids_verify should be allocated")
    call assert_equal_array_char(gene_ids, gene_ids_verify, 128, size(gene_ids), "Gene IDs should match")

    call assert_true(allocated(kallisto_verify), "kallisto_verify should be allocated")
    call assert_equal_array_real(kallisto_expr, kallisto_verify, size(kallisto_expr), 1e-12_real64, "Expression data should match")

    call assert_true(allocated(gene_to_fam_verify), "Gene to family mapping should be allocated")
    call assert_equal_array_int(gene_to_fam, gene_to_fam_verify, size(gene_to_fam), "Gene to family mapping should match")

    call assert_true(allocated(gene_family_ids_verify), "Family IDs should be allocated")
    call assert_equal_array_char(gene_family_ids, gene_family_ids_verify, 128, size(gene_family_ids), "Family IDs should match")

    call assert_true(allocated(family_centroids_verify), "Family centroids should be allocated")
    call assert_equal_array_real(family_centroids, family_centroids_verify, size(family_centroids), 1e-12_real64, "Family centroids should match")

    call assert_true(allocated(shift_vectors_verify), "Shift vectors should be allocated")
    call assert_equal_array_real(shift_vectors, shift_vectors_verify, size(shift_vectors), 1e-12_real64, "Shift vectors should match")
    
    ! Clean up
    deallocate(gene_ids_verify)
    deallocate(kallisto_verify)
    deallocate(gene_to_fam_verify)
    deallocate(gene_family_ids_verify)
    deallocate(family_centroids_verify)
    deallocate(shift_vectors_verify)
    
    ! Test 2: Save only gene_ids and expression
    ! print *, "Test 2: Saving only gene_ids and expression"
    call save_tox_data("test_archive_2_f.zip", ierr, &
                      gene_ids=gene_ids, gene_ids_file="gene_ids_v2.bin", &
                      expression=kallisto_expr, expression_file="kallisto_v2.bin")
    
    call assert_equal_int(ierr, ERR_OK, "Error saving archive")
    
    call read_tox_data("test_archive_2_f.zip", ierr, &
                      gene_ids=gene_ids_verify, &
                      expression=kallisto_verify)
    
    call assert_equal_int(ierr, ERR_OK, "Error reading archive")
    
    ! Verify the data
    call assert_true(allocated(gene_ids_verify), "Gene IDs should be allocated")
    call assert_equal_array_char(gene_ids, gene_ids_verify, 128, size(gene_ids), "Gene IDs should match")

    call assert_true(allocated(kallisto_verify), "Expression data should be allocated")
    call assert_equal_array_real(kallisto_expr, kallisto_verify, size(kallisto_expr), 1e-12_real64, "Expression data should match")
    
    ! Try to read arrays that weren't saved (should not be allocated)
    call read_tox_data("test_archive_2_f.zip", ierr, &
                      gene_to_family=gene_to_fam_verify, &
                      family_ids=gene_family_ids_verify, &
                      family_centroids=family_centroids_verify, &
                      shift_vectors=shift_vectors_verify)
    
    call assert_equal_int(ierr, ERR_OK, "Error reading archive for missing arrays")

    call assert_true(allocated(gene_ids_verify), "gene_ids_verify should be allocated")
    call assert_true(allocated(kallisto_verify), "kallisto_verify should be allocated")
    
    call assert_false(allocated(gene_to_fam_verify), "gene_to_fam_verify should not be allocated")

    call assert_false(allocated(gene_family_ids_verify), "gene_family_ids_verify should not be allocated")
    
    call assert_false(allocated(family_centroids_verify), "family_centroids_verify should not be allocated")
    
    call assert_false(allocated(shift_vectors_verify), "shift_vectors_verify should not be allocated")
    
    ! Clean up
    deallocate(gene_ids_verify)
    deallocate(kallisto_verify)
    
    ! Test 3: Save only family data
    ! print *, "Test 3: Saving only family data"
    call save_tox_data("test_archive_3_f.zip", ierr, &
                      family_ids=gene_family_ids, family_ids_file="family_ids_v3.bin", &
                      family_centroids=family_centroids, family_centroids_file="family_centroids_v3.bin")
    
    call assert_equal_int(ierr, ERR_OK, "Error saving archive")
    
    call read_tox_data("test_archive_3_f.zip", ierr, &
                      family_ids=gene_family_ids_verify, &
                      family_centroids=family_centroids_verify)
    
    call assert_equal_int(ierr, ERR_OK, "Error reading archive")
    
    call assert_true(allocated(gene_family_ids_verify), "Family IDs should be allocated")
    call assert_equal_array_char(gene_family_ids, gene_family_ids_verify, 128, size(gene_family_ids), "Family IDs should match")

    call assert_true(allocated(family_centroids_verify), "Family centroids should be allocated")
    call assert_equal_array_real(family_centroids, family_centroids_verify, size(family_centroids), 1e-12_real64, "Family centroids should match")
    
    ! Clean up
    deallocate(gene_family_ids_verify)
    deallocate(family_centroids_verify)
    
    ! Test 4: Save empty archive (should work without error)
    ! print *, "Test 4: Saving empty archive"
    call save_tox_data("test_archive_4_f.zip", ierr)
    
    call assert_equal_int(ierr, ERR_OK, "Error saving empty archive")
    
    call read_tox_data("test_archive_4_f.zip", ierr)
    
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Reading empty archive should throw error 201")
    
    ! print *, "Empty archive test passed"
    
    ! Test 5: Try to read non-existent archive
    ! print *, "Test 5: Trying to read non-existent archive"
    call read_tox_data("non_existent.zip", ierr)
    
    call assert_not_equal_int(ierr, ERR_OK, "Reading non-existent archive should fail")

    ! print *, "Reading R archive"
    call read_tox_data("test_archive_1_R.zip", ierr, &
                      gene_ids=gene_ids_verify, &
                      expression=kallisto_verify, &
                      gene_to_family=gene_to_fam_verify, &
                      family_ids=gene_family_ids_verify, &
                      family_centroids=family_centroids_verify, &
                      shift_vectors=shift_vectors_verify)
    if(is_ok(ierr)) then
      call assert_true(allocated(gene_ids_verify), "Gene IDs should be allocated from R archive")
      call assert_equal_array_char(gene_ids, gene_ids_verify, 128, size(gene_ids), "R Gene IDs should match")

      call assert_true(allocated(kallisto_verify), "Expression data should be allocated from R archive")
      call assert_equal_array_real(kallisto_expr, kallisto_verify, size(kallisto_expr), 1e-12_real64, "R Expression data should match")

      call assert_true(allocated(gene_to_fam_verify), "Gene to family mapping should be allocated from R archive")
      call assert_equal_array_int(gene_to_fam, gene_to_fam_verify, size(gene_to_fam), "R Gene to family mapping should match")

      call assert_true(allocated(gene_family_ids_verify), "Family IDs should be allocated from R archive")
      call assert_equal_array_char(gene_family_ids, gene_family_ids_verify, 128, size(gene_family_ids), "R Family IDs should match")

      call assert_true(allocated(family_centroids_verify), "Family centroids should be allocated from R archive")
      call assert_equal_array_real(family_centroids, family_centroids_verify, size(family_centroids), 1e-12_real64, "R Family centroids should match")

      call assert_true(allocated(shift_vectors_verify), "Shift vectors should be allocated from R archive")
      call assert_equal_array_real(shift_vectors, shift_vectors_verify, size(shift_vectors), 1e-12_real64, "R Shift vectors should match")
      deallocate(gene_ids_verify)
      deallocate(kallisto_verify)
      deallocate(gene_to_fam_verify)
      deallocate(gene_family_ids_verify)
      deallocate(family_centroids_verify)
      deallocate(shift_vectors_verify)
    
      else
      write(*,*) 'Error reading R archive: ', ierr 
      call set_ok(ierr)
    end if
    

    call read_tox_data("test_archive_1_R.zip", ierr, &
                      gene_ids=gene_ids_verify, &
                      expression=kallisto_verify)
    if(is_ok(ierr)) then
      call assert_true(allocated(gene_ids_verify), "Gene IDs should be allocated from R archive")
      call assert_equal_array_char(gene_ids, gene_ids_verify, 128, size(gene_ids), "R Gene IDs should match")

      call assert_true(allocated(kallisto_verify), "Expression data should be allocated from R archive")
      call assert_equal_array_real(kallisto_expr, kallisto_verify, size(kallisto_expr), 1e-12_real64, "R Expression data should match")

      call assert_false(allocated(gene_to_fam_verify), "Gene to family mapping should not be allocated from R archive")

      call assert_false(allocated(gene_family_ids_verify), "Family IDs should not be allocated from R archive")

      call assert_false(allocated(family_centroids_verify), "Family centroids should not be allocated from R archive")

      call assert_false(allocated(shift_vectors_verify), "Shift vectors should not be allocated from R archive")

      deallocate(gene_ids_verify)
      deallocate(kallisto_verify)
    else
      write(*,*) 'Error reading R archive: ', ierr 
      call set_ok(ierr)
    end if



    ! print *, "Reading python archive"
    call read_tox_data("test_archive_1_py.zip", ierr, &
                      gene_ids=gene_ids_verify, &
                      expression=kallisto_verify, &
                      gene_to_family=gene_to_fam_verify, &
                      family_ids=gene_family_ids_verify, &
                      family_centroids=family_centroids_verify, &
                      shift_vectors=shift_vectors_verify)
    if(is_ok(ierr)) then
      call assert_true(allocated(gene_ids_verify), "Gene IDs should be allocated from python archive")
      call assert_equal_array_char(gene_ids, gene_ids_verify, 128, size(gene_ids), "python Gene IDs should match")

      call assert_true(allocated(kallisto_verify), "Expression data should be allocated from python archive")
      call assert_equal_array_real(kallisto_expr, kallisto_verify, size(kallisto_expr), 1e-12_real64, "python Expression data should match")

      call assert_true(allocated(gene_to_fam_verify), "Gene to family mapping should be allocated from python archive")
      call assert_equal_array_int(gene_to_fam, gene_to_fam_verify, size(gene_to_fam), "python Gene to family mapping should match")

      call assert_true(allocated(gene_family_ids_verify), "Family IDs should be allocated from python archive")
      call assert_equal_array_char(gene_family_ids, gene_family_ids_verify, 128, size(gene_family_ids), "python Family IDs should match")

      call assert_true(allocated(family_centroids_verify), "Family centroids should be allocated from python archive")
      call assert_equal_array_real(family_centroids, family_centroids_verify, size(family_centroids), 1e-12_real64, "python Family centroids should match")

      call assert_true(allocated(shift_vectors_verify), "Shift vectors should be allocated from python archive")
      call assert_equal_array_real(shift_vectors, shift_vectors_verify, size(shift_vectors), 1e-12_real64, "python Shift vectors should match")
      deallocate(gene_ids_verify)
      deallocate(kallisto_verify)
      deallocate(gene_to_fam_verify)
      deallocate(gene_family_ids_verify)
      deallocate(family_centroids_verify)
      deallocate(shift_vectors_verify)
    else
      write(*,*) 'Error reading Python archive: ', ierr 
      call set_ok(ierr)
    end if



    ! print *, "All archive tests completed successfully!"
  end subroutine test_archive

  subroutine test_manual_archive()
    integer(int32) :: ierr
    integer(int32) :: i, j
    integer(int32) :: dims(5)
    integer(int32) :: ndims, clen
    
    ! Test arrays
    integer(int32) :: int_1d(5)
    integer(int32) :: int_2d(3, 4)
    real(real64) :: real_1d(6)
    real(real64) :: real_2d(2, 3)
    character(len=20) :: char_1d(4)
    
    ! Filenames
    character(len=32) :: keys(5)
    character(len=256) :: filenames(5)
    
    ! For extraction
    character(len=:), allocatable :: extracted_keys(:)
    character(len=:), allocatable :: extracted_filenames(:)
    
    ! Arrays for verification
    integer(int32), allocatable :: read_int_1d(:)
    integer(int32), allocatable :: read_int_2d(:,:)
    real(real64), allocatable :: read_real_1d(:)
    real(real64), allocatable :: read_real_2d(:,:)
    character(len=:), allocatable :: read_char_1d(:)
    
    call set_ok(ierr)
    
    print *, "=== Starting manual archive test ==="
    
    ! Initialize test arrays
    print *, "Initializing test arrays..."
    
    ! 1D integer array
    do i = 1, 5
        int_1d(i) = i * 10
    end do
    
    ! 2D integer array
    do j = 1, 4
        do i = 1, 3
            int_2d(i, j) = i * 100 + j
        end do
    end do
    
    ! 1D real array
    do i = 1, 6
        real_1d(i) = i * 1.5_real64
    end do
    
    ! 2D real array
    do j = 1, 3
        do i = 1, 2
            real_2d(i, j) = i * 2.5_real64 + j * 0.1_real64
        end do
    end do
    
    ! 1D character array
    char_1d(1) = "Hello"
    char_1d(2) = "World"
    char_1d(3) = "Fortran"
    char_1d(4) = "Testing"
    
    ! Serialize arrays to files
    print *, "Serializing arrays to files..."
    
    call serialize_int_1d(int_1d, "test_int_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error serializing int_1d")
    
    call serialize_int_2d(int_2d, "test_int_2d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error serializing int_2d")
    
    call serialize_real_1d(real_1d, "test_real_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error serializing real_1d")
    
    call serialize_real_2d(real_2d, "test_real_2d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error serializing real_2d")
    
    call serialize_char_1d(char_1d, "test_char_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error serializing char_1d")
    
    ! Set up keys and filenames for ZIP archive
    keys(1) = 'integer_1d'
    filenames(1) = 'test_int_1d.bin'
    
    keys(2) = 'integer_2d'
    filenames(2) = 'test_int_2d.bin'
    
    keys(3) = 'real_1d'
    filenames(3) = 'test_real_1d.bin'
    
    keys(4) = 'real_2d'
    filenames(4) = 'test_real_2d.bin'
    
    keys(5) = 'character_1d'
    filenames(5) = 'test_char_1d.bin'
    
    ! Create ZIP archive
    print *, "Creating ZIP archive..."
    call create_zip_archive("test_archive_manual_1.zip", keys, filenames, ierr)
    call assert_equal_int(ierr, ERR_OK, "Error creating ZIP archive")
    
    print *, "ZIP archive created successfully"
    
    ! Delete original files to test extraction
    ! print *, "Deleting original files..."
    ! call delete_file("test_int_1d.bin", ierr)
    ! call delete_file("test_int_2d.bin", ierr)
    ! call delete_file("test_real_1d.bin", ierr)
    ! call delete_file("test_real_2d.bin", ierr)
    ! call delete_file("test_char_1d.bin", ierr)
    
    ! Extract ZIP archive
    print *, "Extracting ZIP archive..."
    call extract_zip_archive("test_archive_manual_1.zip", extracted_keys, extracted_filenames, ierr)
    call assert_equal_int(ierr, ERR_OK, "Error extracting ZIP archive")
    
    !print *, "ZIP archive extracted successfully"
    
    ! Print extracted keys and filenames
    !print *, "Extracted files:"
    do i = 1, size(extracted_keys)
        print *, "  Key: ", trim(extracted_keys(i)), " -> File: ", trim(extracted_filenames(i))
    end do
    
    ! Verify that files exist and can be read
    !print *, "Verifying extracted files..."
    
    ! Read and verify 1D integer array
    call get_array_metadata('test_int_1d.bin', dims, 1, ndims, ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading metadata for 1D integer array")

    allocate(read_int_1d(dims(1)))
    call deserialize_int_1d(read_int_1d, "test_int_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading 1D integer array")
    
    ! Read and verify 2D integer array
    call get_array_metadata('test_int_2d.bin', dims, 2, ndims, ierr)
    allocate(read_int_2d(dims(1), dims(2)))
    call deserialize_int_2d(read_int_2d, "test_int_2d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading 2D integer array")
    
    ! Read and verify 1D real array
    call get_array_metadata('test_real_1d.bin', dims, 1, ndims, ierr)
    allocate(read_real_1d(dims(1)))
    call deserialize_real_1d(read_real_1d, "test_real_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading 1D real array")
    
    ! Read and verify 2D real array
    call get_array_metadata('test_real_2d.bin', dims, 2, ndims, ierr)
    allocate(read_real_2d(dims(1), dims(2)))
    call deserialize_real_2d(read_real_2d, "test_real_2d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading 2D real array")
    
    ! Read and verify 1D character array
    call get_array_metadata('test_char_1d.bin', dims, 1, ndims, ierr, clen)
    allocate(character(len=clen) :: read_char_1d(dims(1)))
    call deserialize_char_1d(read_char_1d, "test_char_1d.bin", ierr)
    call assert_equal_int(ierr, ERR_OK, "Error reading 1D character array")
    
    ! Clean up extracted files
    ! call delete_file("test_inc
    
    ! Clean up extracted files
    ! call delete_file("test_int_1d.bin", ierr)
    ! call delete_file("test_int_2d.bin", ierr)
    ! callt_1d.bin", ierr)
    ! call delete_file("test_int_2d.bin", ierr)
    ! call delete_file("test_real_1d.bin", ierr)
    ! call delete_file("test_real_2d.bin", ierr)
    ! call delete_file("test_char_1d.bin", ierr)
    ! call delete_file("test_manual_archive.zip", ierr)
    ! call delete_file("manifest.txt", ierr)
  end subroutine test_manual_archive

  subroutine test_hashing()
    use f42_xxh3_hashmap
    type(hashmap_type) :: test_hashmap
    type(hashset_type) :: test_hashset
    character(len=6), allocatable :: keys(:)
    integer(int32), allocatable :: values(:)
    integer(int32) :: value, i, ierr
    logical :: in_hashset

    allocate(keys(5))
    allocate(values(5))
    call set_ok(ierr)

    keys = [ &
      'First ', &
      'Second', &
      'Foo   ', &
      'Bar   ', &
      'First ' &
    ]

    values = [1,2,3,4,5]

    call hashmap_create(test_hashmap)

    do i = 1, 4
      call hashmap_put(test_hashmap, keys(i), values(i))
    end do 
    value = hashmap_get(test_hashmap, 'First')
    call assert_equal_int(value, 1, 'Wrong return value for key First')

    value = hashmap_get(test_hashmap, 'Second')
    call assert_equal_int(value, 2, 'Wrong return value for key Second')

    value = hashmap_get(test_hashmap, 'Foo')
    call assert_equal_int(value, 3, 'Wrong return value for key Foo')

    value = hashmap_get(test_hashmap, 'Bar')
    call assert_equal_int(value, 4, 'Wrong return value for key Bar')

    call hashmap_put(test_hashmap, 'First', 5)
    value = hashmap_get(test_hashmap, 'First')
    call assert_equal_int(value, 5, 'Wrong return for updated key First')

    call hashset_create(test_hashset)
    do i = 1, 5
      call hashset_put(test_hashset, keys(i), ierr)
    end do
    call assert_equal_int(ierr, ERR_INVALID_INPUT, 'Inserting duplicate key should return error')
    call set_ok(ierr)

    do i = 1, 5
      in_hashset = is_in_hashset(test_hashset, keys(i))
      call assert_true(in_hashset, 'Key should be in hashset')
    end do

    in_hashset = is_in_hashset(test_hashset, 'Dummy')
    call assert_false(in_hashset, 'Dummy should not be in hashset')

    call hashset_destroy(test_hashset)
    call hashmap_destroy(test_hashmap)
  end subroutine test_hashing

  !> Test duplicate gene IDs detection
  subroutine test_validate_duplicate_gene_ids()
      character(len=20), allocatable :: test_gene_ids(:)
      integer(int32) :: ierr
      
      allocate(test_gene_ids(3))
      test_gene_ids = ['NP_001000001.1', 'NP_001000002.1', 'NP_001000001.1']  ! Duplicate
      
      call validate_string_array_uniqueness(test_gene_ids, ierr)
      call assert_equal_int(ierr, ERR_INVALID_INPUT, "Should detect duplicate gene IDs")
      
      deallocate(test_gene_ids)
  end subroutine test_validate_duplicate_gene_ids

  !> Test NaN/Inf detection in expression data
  !> Test NaN detection in validation (bypassing file reading checks)
  subroutine test_validate_nan_direct()
    use, intrinsic :: ieee_arithmetic
    real(real64), allocatable :: test_expr(:,:)
    integer(int32) :: ierr
    
    allocate(test_expr(2, 3))
    
    ! Fill with normal data first
    test_expr = reshape([1.0_real64, 2.0_real64, 3.0_real64, &
                        4.0_real64, 5.0_real64, 6.0_real64], [2, 3])
    
    ! Inject NaN directly into the array (simulating computational error)
    test_expr(2, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
    
    ! Test the validation routine directly
    call validate_expression_data(test_expr, .false., ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Validation should detect computational NaN")
    
    deallocate(test_expr)
  end subroutine test_validate_nan_direct

  !> Test Inf detection in validation  
  subroutine test_validate_inf_direct()
    use, intrinsic :: ieee_arithmetic
    real(real64), allocatable :: test_expr(:,:)
    integer(int32) :: ierr
    
    allocate(test_expr(2, 3))
    
    test_expr = reshape([1.0_real64, 2.0_real64, 3.0_real64, &
                        4.0_real64, 5.0_real64, 6.0_real64], [2, 3])
    
    ! Inject Inf directly
    test_expr(1, 3) = ieee_value(1.0_real64, ieee_positive_inf)
    
    call validate_expression_data(test_expr, .false., ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Validation should detect computational Inf")
    
    deallocate(test_expr)
  end subroutine test_validate_inf_direct

  !> Test negative expression values
  subroutine test_validate_negative_expression()
      real(real64), allocatable :: test_expr(:,:)
      integer(int32) :: ierr
      
      allocate(test_expr(2, 2))
      test_expr = reshape([1.0_real64, -2.0_real64, 3.0_real64, 4.0_real64], [2,2])  ! Negative value
      
      call validate_expression_data(test_expr, .true., ierr)
      call assert_equal_int(ierr, ERR_INVALID_INPUT, "Should detect negative expression values")
      
      deallocate(test_expr)
  end subroutine test_validate_negative_expression

  !> Test invalid gene-to-family mapping
  subroutine test_validate_invalid_gene_to_fam()
      integer(int32), allocatable :: test_gene_to_fam(:)
      integer(int32) :: ierr, n_families
      
      allocate(test_gene_to_fam(3))
      test_gene_to_fam = [1, 2, 5]  ! Family 5 doesn't exist
      n_families = 3  ! Only families 1,2,3 exist
      
      call validate_gene_to_family_mapping(test_gene_to_fam, n_families, ierr)
      call assert_equal_int(ierr, ERR_INVALID_INPUT, "Should detect invalid family indices")
      
      deallocate(test_gene_to_fam)
  end subroutine test_validate_invalid_gene_to_fam

  !> Test empty strings validation
  subroutine test_validate_empty_strings()
      character(len=20), allocatable :: test_strings(:)
      integer(int32) :: ierr
      
      allocate(test_strings(3))
      test_strings = ['NP_001000001.1', '              ', 'NP_001000003.1']  ! Empty string
      
      call validate_empty_strings(test_strings, "test_array", ierr)
      call assert_equal_int(ierr, ERR_INVALID_INPUT, "Should detect empty strings")
      
      deallocate(test_strings)
  end subroutine test_validate_empty_strings

  !> Test dimension mismatch validation
  subroutine test_validate_dimension_mismatch()
      integer(int32) :: ierr
      integer(int32) :: n_genes, n_families, n_samples
      character(len=20), allocatable :: gene_ids(:), gene_family_ids(:)
      integer(int32), allocatable :: gene_to_fam(:)
      real(real64), allocatable :: expression_vectors(:,:), family_centroids(:,:), shift_vectors(:,:)
      
      ! Setup with intentional mismatch
      n_genes = 2
      n_families = 1
      n_samples = 2
      
      allocate(gene_ids(3))  ! Wrong size - should be n_genes=2
      allocate(gene_family_ids(n_families))
      allocate(gene_to_fam(n_genes))
      allocate(expression_vectors(n_samples, n_genes))
      allocate(family_centroids(n_samples, n_families))
      allocate(shift_vectors(2*n_samples, n_genes))
      
      call validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, &
                                gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
      call assert_equal_int(ierr, ERR_SIZE_MISMATCH, "Should detect dimension mismatch")
      
      deallocate(gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors)
  end subroutine test_validate_dimension_mismatch
end module mod_test_tox_data