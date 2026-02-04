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
    all_tests(2) = test_case("test_rdi_basic", test_rdi_basic)
    all_tests(3) = test_case("test_identify_outliers_basic", test_identify_outliers_basic)
    all_tests(4) = test_case("test_detect_outliers_basic", test_detect_outliers_basic)
    all_tests(5) = test_case("test_detect_outliers_invalid_indices", test_detect_outliers_invalid_indices)
    all_tests(6) = test_case("test_all_outliers", test_all_outliers)
    all_tests(7) = test_case("test_no_outliers", test_no_outliers)
    all_tests(8) = test_case("test_percentile_boundaries", test_percentile_boundaries)
    all_tests(9) = test_case("test_identical_rdi_at_threshold", test_identical_rdi_at_threshold)
    all_tests(10) = test_case("test_all_negative_rdi", test_all_negative_rdi)
    all_tests(11) = test_case("test_single_gene_family", test_single_gene_family)
    all_tests(12) = test_case("test_invalid_indices", test_invalid_indices)
    all_tests(13) = test_case("test_all_zero_distances", test_all_zero_distances)
    all_tests(14) = test_case("test_all_genes_one_family", test_all_genes_one_family)
    all_tests(15) = test_case("test_negative_distances", test_negative_distances)
    all_tests(16) = test_case("test_invalid_family_indices", test_invalid_family_indices)
    all_tests(17) = test_case("test_loess_global_positive", test_loess_global_positive)
    all_tests(18) = test_case("test_single_gene_family_scaling", test_single_gene_family_scaling)
    all_tests(19) = test_case("test_all_zero_distances_scaling", test_all_zero_distances_scaling)

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
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err
    implicit none

    integer(int32), parameter :: n_families = 6_int32
    integer(int32), parameter :: genes_per_fam = 4_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: ierr
    real(real64), parameter :: offsets(genes_per_fam) = [0.1_real64, 0.4_real64, 0.9_real64, 1.6_real64]


    integer(int32) :: f, j, idx
    integer(int32) :: n_valid

    ! Build test data: 4 genes per family, with varying distances to ensure stddev > 0
    idx = 0
    do f = 1, n_families
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = f
        distances(idx) = real(f, real64) * 2.0_real64 + offsets(j)
      end do
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 (test_scaling_basic)')
    call assert_true(all(dscale >= 0.0_real64), 'Scaling should be non-negative (test_scaling_basic)')

    ! verify that all families were considered valid (since each has 4 genes)
    n_valid = 0
    do j = 1, n_families
      if (indices_used(j) >= 1_int32 .and. indices_used(j) <= n_families) n_valid = n_valid + 1
    end do
    call assert_equal_int(n_valid, n_families, 'All families should be valid (test_scaling_basic)')

    ! basic sanity that LOESS x/y got populated for valid families
    call assert_true(all(loess_x > 0.0_real64), 'loess_x should be > 0 for valid families')
    call assert_true(all(loess_y > 0.0_real64), 'loess_y should be > 0 for valid families')
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
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 6_int32
    integer(int32), parameter :: genes_per_fam = 6_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical        :: is_outlier(n_genes)
    integer(int32) :: ierr

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)

    integer(int32) :: f, j, idx
    real(real64), parameter :: offsets(genes_per_fam) = [0.10_real64, 0.30_real64, 0.50_real64, 0.80_real64, 1.10_real64, 20.0_real64]
    ! last entry is a deliberate "big" outlier within each family

    idx = 0
    do f = 1, n_families
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = f
        distances(idx) = real(f, real64) * 2.0_real64 + offsets(j)
      end do
    end do

    call detect_outliers( &
        n_genes, n_families, distances, gene_to_fam, &
        work_array, perm, stack_left, stack_right, &
        is_outlier, loess_x, loess_y, loess_n, ierr, &
        60.0_real64)

    call assert_equal_int(ierr, 0, "ierr should be 0 (test_detect_outliers_basic)")
    call assert_true(any(is_outlier), "At least one outlier detected")
  end subroutine test_detect_outliers_basic



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


  !> Test: Multiple genes with identical RDI at threshold.
  !> This test checks that if multiple genes have identical RDI at the threshold, they are all identified as outliers.
  subroutine test_identical_rdi_at_threshold()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [1, 2, 2, 3] already sorted and no negatives
    sorted_rdi = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 50.0_real64)
    call assert_true(is_outlier(3), "Identical RDI at threshold is outlier")
    call assert_true(is_outlier(2), "Identical RDI at threshold is outlier")
  end subroutine test_identical_rdi_at_threshold


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
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, ierr)
    call assert_equal_int(ierr, 201, 'detect_outliers propagates error_code 201 for invalid indices')
  end subroutine test_detect_outliers_invalid_indices


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
    ! rdi = [-0.1, -0.2, -0.3, -0.4, -0.5] -> sorted_rdi = [0,0,0,0,0] (all negatives)
    sorted_rdi = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 80.0_real64)
    call assert_true(.not. any(is_outlier), "All-negative RDI should not be outliers")
  end subroutine test_all_negative_rdi

  subroutine test_invalid_indices()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 10_int32
    integer(int32), parameter :: n_genes    = 50_int32

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: ierr

    integer(int32) :: i

    ! Fill data
    do i = 1, n_genes
      distances(i) = real(i, real64)
      gene_to_fam(i) = 1_int32 + mod(i-1, n_families)  ! valid 1..n_families
    end do

    ! Inject invalid indices
    gene_to_fam(7)  = 0_int32
    gene_to_fam(23) = n_families + 3_int32

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 201, 'Error code 201 for invalid family indices')
    call assert_true(all(dscale == -1.0_real64), 'dscale must be -1.0 on error')
  end subroutine test_invalid_indices

  !> Test: Single gene and family case.
  !> This test checks that for a single gene in a single family, scaling and RDI are both 0.0, and the gene is not an outlier.
  subroutine test_single_gene_family()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes = 1_int32, n_families = 1_int32

    real(real64)   :: distances(n_genes)   = [0.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1_int32]

    real(real64)   :: dscale(n_families)
    real(real64)   :: rdi(n_genes)
    real(real64)   :: sorted_rdi(n_genes)
    real(real64)   :: work_array(n_genes)

    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical        :: is_outlier(n_genes)
    integer(int32) :: ierr
    integer(int32) :: i

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)  ! NOTE: used as workspace by scaling alloc in your detect_outliers

    ! 1) scaling (alloc wrapper)
    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, loess_n, ierr)
    call assert_equal_int(ierr, 0, "Scaling ierr should be 0 (test_single_gene_family)")

    ! 2) RDI
    perm = [(i, i=1,n_genes)]
    stack_left  = 0_int32
    stack_right = 0_int32

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! 3) Outlier detection (uses scaling alloc inside detect_outliers too, but that's fine)
    perm = [(i, i=1,n_genes)]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                        work_array, perm, stack_left, stack_right, &
                        is_outlier, loess_x, loess_y, loess_n, ierr, &
                        95.0_real64)
    call assert_equal_int(ierr, 0, "detect_outliers ierr should be 0 (test_single_gene_family)")

    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, "Single gene scaling (should be 0.0)")
    call assert_equal_real(rdi(1),    0.0_real64, 1e-12_real64, "Single gene RDI")
    call assert_false(is_outlier(1), "Single gene not outlier")
  end subroutine test_single_gene_family

  !> Test: All distances are zero.
  !> This test checks that if all distances are zero, scaling is 0.0 and all RDI values are 0.0.
  subroutine test_all_zero_distances()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 6_int32
    integer(int32), parameter :: genes_per_fam = 5_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: rdi(n_genes), sorted_rdi(n_genes)

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)     ! workspace for scaling alloc (indices_used)

    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    integer(int32) :: ierr
    integer(int32) :: i, f, j, idx

    ! All distances = 0
    distances = 0.0_real64

    ! Map genes to families evenly
    idx = 0
    do f = 1, n_families
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = f
      end do
    end do

    ! Scaling (alloc wrapper)
    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, loess_n, ierr)
    call assert_equal_int(ierr, 0, "Scaling ierr should be 0 (test_all_zero_distances)")

    ! Prepare sort work arrays for compute_rdi
    perm = [(i, i=1,n_genes)]
    stack_left  = 0_int32
    stack_right = 0_int32

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! Assertions
    call assert_true(all(dscale == 0.0_real64), "All zero distances => all dscale == 0.0")
    call assert_true(all(rdi    == 0.0_real64), "All zero distances => all RDI == 0.0")
  end subroutine test_all_zero_distances

  subroutine test_all_genes_one_family()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 1_int32
    integer(int32), parameter :: n_genes    = 40_int32

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    integer(int32) :: ierr

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)

    integer(int32) :: i

    ! All genes in family 1
    gene_to_fam = 1_int32

    ! Give varied distances (non-trivial family), but LOESS-global is undefined with only 1 family
    do i = 1, n_genes
      distances(i) = 1.0_real64 + 0.1_real64 * real(i-1, real64)
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 0, "ierr should be 0 (test_all_genes_one_family)")

    ! With only one valid family (n_valid==1), compute_family_scaling returns early and leaves dscale = 0.
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, &
        "One family only => n_valid==1 => dscale should remain 0.0")
  end subroutine test_all_genes_one_family

  subroutine test_negative_distances()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 12_int32
    integer(int32), parameter :: n_genes    = 12_int32  ! one gene per family => all singleton

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: rdi(n_genes), sorted_rdi(n_genes)

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)  ! workspace for scaling alloc

    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    integer(int32) :: ierr
    integer(int32) :: i

    ! Start with all zeros
    distances = 0.0_real64
    do i = 1, n_genes
      gene_to_fam(i) = i   ! 1 gene per family => singleton families
    end do

    ! Scaling: with all singleton families, n_valid stays 0 => dscale remains 0 for all families
    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)
    call assert_equal_int(ierr, 0, "Scaling ierr should be 0 (test_negative_nan_distances)")
    call assert_true(all(dscale == 0.0_real64), "All singleton families => dscale should remain 0.0")

    ! Inject a negative distance AFTER scaling (we are testing compute_rdi behavior)
    distances(1) = -1.0_real64

    perm = [(i, i=1,n_genes)]
    stack_left  = 0_int32
    stack_right = 0_int32

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    call assert_equal_real(rdi(1), 0.0_real64, 1e-12_real64, &
        "Negative distance with zero scaling gives RDI=0.0 (not outlier)")

  end subroutine test_negative_distances

  !> Exhaustive: invalid family indices (should error 201).
  !> This test checks that if gene_to_fam contains invalid indices, error_code is 201 and dscale falls back to -1.0 for all families.
  subroutine test_invalid_family_indices()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 12_int32
    integer(int32), parameter :: n_genes    = 60_int32

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    integer(int32) :: ierr

    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)

    integer(int32) :: i

    ! Build valid data first
    do i = 1, n_genes
      distances(i)   = real(i, real64)
      gene_to_fam(i) = 1_int32 + mod(i-1, n_families)   ! valid 1..n_families
    end do

    ! Inject a couple invalid indices
    gene_to_fam(10) = 0_int32
    gene_to_fam(37) = n_families + 5_int32

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 201, 'Error code 201 for invalid family indices')
    call assert_true(all(dscale == -1.0_real64), 'dscale fallback to -1.0 on invalid indices')
  end subroutine test_invalid_family_indices

  subroutine test_single_gene_family_scaling()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 12_int32
    integer(int32), parameter :: n_genes    = 12_int32  ! one gene per family

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: ierr
    integer(int32) :: i

    do i = 1, n_genes
      distances(i)   = real(i, real64) * 1.1_real64
      gene_to_fam(i) = i  ! singleton families
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 for single-gene families')
    call assert_true(all(dscale == 0.0_real64), 'All singleton families scaling 0.0')
  end subroutine test_single_gene_family_scaling

  subroutine test_all_zero_distances_scaling()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 6_int32
    integer(int32), parameter :: genes_per_fam = 6_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: ierr

    integer(int32) :: f, j, idx

    distances = 0.0_real64

    idx = 0
    do f = 1, n_families
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = f
      end do
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 for all-zero distances')
    call assert_true(all(dscale == 0.0_real64), 'All-zero distances scaling 0.0')
  end subroutine test_all_zero_distances_scaling

  subroutine test_loess_global_positive()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 8_int32
    integer(int32), parameter :: genes_per_fam = 6_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)   :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)

    real(real64)   :: dscale(n_families)
    real(real64)   :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: ierr

    integer(int32) :: f, j, idx
    real(real64), parameter :: offsets(genes_per_fam) = [0.1_real64, 0.4_real64, 0.9_real64, 1.5_real64, 2.2_real64, 3.0_real64]

    idx = 0
    do f = 1, n_families
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = f
        distances(idx)   = real(f, real64) * 2.5_real64 + offsets(j)
      end do
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 for LOESS global scaling')
    call assert_true(any(dscale > 0.0_real64), 'At least one scaling value should be > 0 for non-degenerate input')
  end subroutine test_loess_global_positive


end module mod_test_get_outliers


