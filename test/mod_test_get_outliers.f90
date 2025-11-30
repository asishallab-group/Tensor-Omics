!> @brief Unit test suite for tox_get_outliers routines.
module mod_test_get_outliers
  use asserts
  use tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> @brief Get array of all available tests (as subroutine, not function).
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(19)
    all_tests(1) = test_case("test_scaling_basic", test_scaling_basic)
    all_tests(2) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(3) = test_case("test_rdi_basic", test_rdi_basic)
    all_tests(4) = test_case("test_identify_outliers_basic", test_identify_outliers_basic)
    all_tests(5) = test_case("test_detect_outliers_basic", test_detect_outliers_basic)
    all_tests(6) = test_case("test_all_zero_distances_scaling", test_all_zero_distances_scaling)
    all_tests(7) = test_case("test_detect_outliers_invalid_indices", test_detect_outliers_invalid_indices)
    all_tests(8) = test_case("test_invalid_indices", test_invalid_indices)
    all_tests(9) = test_case("test_all_outliers", test_all_outliers)
    all_tests(10) = test_case("test_no_outliers", test_no_outliers)
    all_tests(11) = test_case("test_single_gene_family", test_single_gene_family)
    all_tests(12) = test_case("test_all_zero_distances", test_all_zero_distances)
    all_tests(13) = test_case("test_all_genes_one_family", test_all_genes_one_family)
    all_tests(14) = test_case("test_percentile_boundaries", test_percentile_boundaries)
    all_tests(15) = test_case("test_negative_nan_distances", test_negative_nan_distances)
    all_tests(16) = test_case("test_identical_rdi_at_threshold", test_identical_rdi_at_threshold)
    all_tests(17) = test_case("test_invalid_family_indices", test_invalid_family_indices)
    all_tests(18) = test_case("test_single_gene_family_scaling", test_single_gene_family_scaling)
    all_tests(19) = test_case("test_all_negative_rdi", test_all_negative_rdi)

    
  end subroutine get_all_tests

  !> @brief Run specific get_outliers tests by name.
  subroutine run_named_tests_get_outliers(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(25)
    integer(int32) :: i, j
    logical :: found

    call get_all_tests(all_tests)

    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_get_outliers

  !> @brief Run all tox_get_outliers tests.
  subroutine run_all_tests_get_outliers()
    type(test_case) :: all_tests(19)
    integer(int32) :: i
    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tox_get_outliers tests passed successfully."
  end subroutine run_all_tests_get_outliers

  !> Test: Basic scaling with orthologs and non-orthologs.
  !> This test checks that the scaling for each family is set to the corresponding value using loess.
  subroutine test_scaling_basic()
    integer(int32), parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,2,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 0, 'Error code 0 (test_scaling_basic)')
    call assert_true(all(dscale > 0.0_real64), 'LOESS scaling should be positive for nontrivial families')
  end subroutine test_scaling_basic

  !> Test: Basic RDI calculation.
  !> This test checks that the RDI (relative distance index) is computed correctly for a simple case with two families.
  subroutine test_rdi_basic()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [2.0_real64, 4.0_real64, 6.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,2,2]
    real(real64) :: dscale(n_families) = [2.0_real64, 4.0_real64]
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    integer(int32) :: i
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_equal_real(rdi(1), 1.0_real64, 1e-12_real64, "RDI gene 1")
    call assert_equal_real(rdi(2), 1.0_real64, 1e-12_real64, "RDI gene 2")
    call assert_equal_real(rdi(3), 1.5_real64, 1e-12_real64, "RDI gene 3")
  end subroutine test_rdi_basic

  !> Test: Outlier detection in a simple RDI vector.
  !> This test checks that the identify_outliers routine correctly identifies the highest and second highest RDI values as outliers at the 75th percentile.
  subroutine test_identify_outliers_basic()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    ! rdi = [0.1, 0.2, 0.3, 0.4] is already sorted and there are no negatives
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 75.0_real64)
    call assert_true(is_outlier(4), "Highest RDI is outlier")
    call assert_true(is_outlier(3), "Second highest RDI is outlier")
    call assert_false(is_outlier(1), "Lowest RDI is not outlier")
  end subroutine test_identify_outliers_basic

  !> Test: Full outlier detection workflow.
  !> This test checks the complete workflow of outlier detection, including scaling, RDI calculation, and outlier identification, for a small example with two families.
  subroutine test_detect_outliers_basic()
    integer(int32), parameter :: n_genes=5, n_families=2
    real(real64) :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    distances = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    gene_to_fam = [1,1,2,2,2]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, ierr, &
      60.0_real64)
    call assert_true(any(is_outlier), "At least one outlier detected")
  end subroutine test_detect_outliers_basic

  !> Test: Invalid family indices.
  !> This test checks that if gene_to_fam contains invalid indices (e.g., 0 or out of range), the error code is set to 201 and dscale falls back to -1.0 for all families.
  subroutine test_invalid_indices()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,3,0] ! 3 and 0 invalid
    real(real64) :: dscale(n_families)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 201, 'Error code 201 for invalid family indices')
    call assert_true(all(dscale == -1.0_real64), 'dscale must be -1.0 on error')
  end subroutine test_invalid_indices

  !> Test: All genes are outliers (percentile 0).
  !> This test checks that if the outlier percentile is set to 0, all genes are identified as outliers.
  subroutine test_all_outliers()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [10, 11, 12, 13] already sorted and no negatives
    sorted_rdi = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 0.0_real64)
    call assert_true(all(is_outlier), "All are outliers at 0 percentile")
  end subroutine test_all_outliers

  !> Test: No genes are outliers (percentile 100).
  !> This test checks that if the outlier percentile is set to 100, no genes except possibly the highest are identified as outliers.
  subroutine test_no_outliers()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [0.1, 0.2, 0.3, 0.4] already sorted and no negatives
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 100.0_real64)
    call assert_true(.not. any(is_outlier(1:3)), "No outliers below 100 percentile")
  end subroutine test_no_outliers

  !> Test: Single gene and family case.
  !> This test checks that for a single gene in a single family, scaling and RDI are both 0.0, and the gene is not an outlier.
  subroutine test_single_gene_family()
    integer(int32), parameter :: n_genes=1, n_families=1
    real(real64) :: distances(1) = [0.0_real64]
    integer(int32) :: gene_to_fam(1) = [1]
    real(real64) :: dscale(1), rdi(1), sorted_rdi(1), work_array(1)
    integer(int32) :: perm(1), stack_left(1), stack_right(1)
    logical :: is_outlier(1)
    integer(int32) :: i, ierr
    real(real64) :: loess_x(1), loess_y(1)
    integer(int32) :: loess_n(1)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, loess_n, perm, stack_left, stack_right, family_distances, ierr)
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    perm = [(i, i=1,n_genes)]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, ierr, &
      95.0_real64)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, &
         "Single gene scaling (should be 0.0)")
    call assert_equal_real(rdi(1), 0.0_real64, 1e-12_real64, &
         "Single gene RDI")
    call assert_false(is_outlier(1), "Single gene not outlier")
  end subroutine test_single_gene_family

  !> Test: All distances are zero.
  !> This test checks that if all distances are zero, scaling is 0.0 and all RDI values are 0.0.
  subroutine test_all_zero_distances()
    integer(int32), parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes) = [0.0_real64, 0.0_real64, 0.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(1), rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: ierr
    real(real64) :: max_distance_bw_orths(1) = [0.0_real64]
    real(real64) :: loess_x(1), loess_y(1), loess_n(1)
    integer(int32) :: indices_used(1)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: i
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, &
         "All zero distances scaling (should be 0.0)")
    call assert_true(all(rdi == 0.0_real64), "All RDI zero for zero distances")
  end subroutine test_all_zero_distances

  !> Test: All genes in a single family.
  !> This test checks that if all genes belong to a single family, scaling is set to the value in max_distance_bw_orths if orthologs are present.
  subroutine test_all_genes_one_family()
    integer(int32), parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1,1]
    real(real64) :: dscale(1)
    integer(int32) :: ierr
    real(real64) :: loess_x(1), loess_y(1), loess_n(1)
    integer(int32) :: indices_used(1)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_true(dscale(1) > 0.0_real64, "All genes one family scaling (LOESS, should be positive)")
  end subroutine test_all_genes_one_family

  !> Test: Percentile at boundaries (0 and 100).
  !> This test checks that at percentile 0, all genes are outliers, and at percentile 100, no genes except possibly the highest are outliers.
  subroutine test_percentile_boundaries()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [1, 2, 3, 4] already sorted and no negatives
    sorted_rdi = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 0.0_real64)
    call assert_true(all(is_outlier), "All outliers at 0 percentile")
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 100.0_real64)
    call assert_true(.not. any(is_outlier(1:3)), "No outliers below 100 percentile")
  end subroutine test_percentile_boundaries

  !> Test: Negative distances only.
  !> This test checks that negative distances with zero scaling result in RDI=0.0 and are not outliers. NaN test is skipped for portability.
  subroutine test_negative_nan_distances()
    integer(int32), parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1]
    real(real64) :: dscale(1), rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: i
    distances = [0.0_real64, 0.0_real64, 0.0_real64]
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    distances(1) = -1.0_real64
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_true(rdi(1) == 0.0_real64, "Negative distance with zero scaling gives RDI=0.0 (not outlier)")
    ! NaN test skipped for portability
  end subroutine test_negative_nan_distances

  !> Test: Multiple genes with identical RDI at threshold.
  !> This test checks that if multiple genes have identical RDI at the threshold, they are all identified as outliers.
  subroutine test_identical_rdi_at_threshold()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    ! rdi = [1, 2, 2, 3] already sorted and no negatives
    sorted_rdi = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 50.0_real64)
    call assert_true(is_outlier(3), "Identical RDI at threshold is outlier")
    call assert_true(is_outlier(2), "Identical RDI at threshold is outlier")
  end subroutine test_identical_rdi_at_threshold

  !> Exhaustive: invalid family indices (should error 201).
  !> This test checks that if gene_to_fam contains invalid indices, error_code is 201 and dscale falls back to -1.0 for all families.
  subroutine test_invalid_family_indices()
    integer(int32), parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,2,3,0]
    real(real64) :: dscale(n_families)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 201, 'Error code 201 for invalid family indices')
    call assert_true(all(dscale == -1.0_real64), 'dscale fallback to -1.0 on invalid indices')
  end subroutine test_invalid_family_indices

  !> Exhaustive: single-gene family (should scale to 0.0).
  !> This test ensures that for families with only one gene, the scaling is set to 0.0 and error code is 0.
  subroutine test_single_gene_family_scaling()
    integer(int32), parameter :: n_genes=2, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 0, 'Error code 0 for single-gene families')
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, 'Single-gene family scaling 0.0')
    call assert_equal_real(dscale(2), 0.0_real64, 1e-12_real64, 'Single-gene family scaling 0.0')
  end subroutine test_single_gene_family_scaling

  !> Exhaustive: all-zero distances (should scale to 0.0).
  !> This test checks that if all distances in a family are zero, the scaling is 0.0 and error code is 0.
  subroutine test_all_zero_distances_scaling()
    integer(int32), parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes) = [0.0_real64, 0.0_real64, 0.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1]
    real(real64) :: dscale(1)
    integer(int32) :: ierr
    real(real64) :: loess_x(1), loess_y(1), loess_n(1)
    integer(int32) :: indices_used(1)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 0, 'Error code 0 for all-zero distances')
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, 'All-zero distances scaling 0.0')
  end subroutine test_all_zero_distances_scaling

  !> Exhaustive: detect_outliers error propagation (invalid indices).
  !> This test checks that detect_outliers propagates error_code 201 if gene_to_fam contains invalid indices.
  subroutine test_detect_outliers_invalid_indices()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,3,0]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: ierr
    logical :: is_ortholog(n_genes) = [.false., .false., .false.]
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, ierr)
    call assert_equal_int(ierr, 201, 'detect_outliers propagates error_code 201 for invalid indices')
  end subroutine test_detect_outliers_invalid_indices

  !> Exhaustive: LOESS fallback for scaling (workspace arrays required).
  subroutine test_loess_fallback()
    integer(int32), parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,2,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families),loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, ierr)
    call assert_equal_int(ierr, 0, 'Error code 0 for LOESS fallback')
    call assert_true(any(dscale > 0.0_real64), 'Scaling must be positive if LOESS fallback is successful')
  end subroutine test_loess_fallback

  !> Test: identify_outliers uses default percentile (95) if not provided.
  subroutine test_identify_outliers_default_percentile()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [0.1, 0.2, 0.3, 0.4] is already sorted and there are no negatives
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold)
    ! At 95th percentile, only the highest value should be outlier
    call assert_true(is_outlier(4), "Highest RDI is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")
  end subroutine test_identify_outliers_default_percentile

  !> Test: detect_outliers uses default percentile (95) if not provided.
  subroutine test_detect_outliers_default_percentile()
    integer(int32), parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1,1]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, ierr)
    ! At 95th percentile, only the highest value should be outlier
    call assert_true(is_outlier(4), "Highest distance is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")
    call assert_equal_int(ierr, 0, "No error for default percentile")
  end subroutine test_detect_outliers_default_percentile

  !> Test: All RDI negative (should be ignored for outlier detection).
  !> This test checks that if all RDI values are negative, none are marked as outliers.
  subroutine test_all_negative_rdi()
    integer(int32), parameter :: n_genes=5
    real(real64) :: rdi(n_genes) = [-0.1_real64, -0.2_real64, -0.3_real64, -0.4_real64, -0.5_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    ! rdi = [-0.1, -0.2, -0.3, -0.4, -0.5] -> sorted_rdi = [0,0,0,0,0] (all negatives)
    sorted_rdi = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 80.0_real64)
    call assert_true(.not. any(is_outlier), "All-negative RDI should not be outliers")
  end subroutine test_all_negative_rdi


end module mod_test_get_outliers


