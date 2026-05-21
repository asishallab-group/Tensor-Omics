!> @brief Unit test suite for tox_get_outliers routines.
module mod_test_get_outliers
  use asserts
  use tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use test_suite
  implicit none
  public

contains

  !> @brief Get array of all available tests (as subroutine, not function).
  function get_all_tests_get_outliers() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(24))

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
    all_tests(20) = test_case("test_outliers_zero_variance", test_outliers_zero_variance)
    all_tests(21) = test_case("test_outliers_orphan_genes", test_outliers_orphan_genes)
    all_tests(22) = test_case("test_outliers_extreme_detection", test_outliers_extreme_detection)
    all_tests(23) = test_case("test_outliers_nan_handling", test_outliers_nan_handling)
    all_tests(24) = test_case("test_scaling_performance_benchmark",test_scaling_performance_benchmark)

  end function get_all_tests_get_outliers

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

  !> Benchmark:9k families and 20k genes
  subroutine test_scaling_performance_benchmark()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err
    implicit none

    integer(int32), parameter :: n_families = 20000_int32
    integer(int32), parameter :: n_genes = 200000_int32

    real(real64), allocatable :: distances(:)
    integer(int32), allocatable :: gene_to_fam(:)
    real(real64), allocatable :: dscale(:)
    real(real64), allocatable :: loess_x(:), loess_y(:)
    integer(int32), allocatable :: indices_used(:)
    
    integer(int32) :: ierr, i, f_idx
    real(real64)    :: t1, t2
    real(real64)    :: dummy_val

    allocate(distances(n_genes), gene_to_fam(n_genes))
    allocate(dscale(n_families), loess_x(n_families), loess_y(n_families))
    allocate(indices_used(n_families))

    ! Distribute 20k genes in 9k families
    do i = 1, n_genes
       ! make sure all families have genes
       f_idx = mod(i, n_families) + 1
       gene_to_fam(i) = f_idx
       ! random distances to make sure we have stddev > 0
       call random_number(dummy_val)
       distances(i) = real(f_idx, real64) * 0.5_real64 + dummy_val
    end do

    print *, "--- Starting Benchmark ---"
    print *, "Genes: ", n_genes
    print *, "Families: ", n_families

    ! 2. Timing subroutine
    call cpu_time(t1)
    
    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)
    
    call cpu_time(t2)

    ! 3. Results
    if (ierr == 0) then
       print *, "Total time: ", (t2 - t1), " seconds."
       print *, "First 5 dscales: ", dscale(1:5)
    else
       print *, "Error detected during scaling: ", ierr
    end if

    ! Clear memory
    deallocate(distances, gene_to_fam, dscale, loess_x, loess_y, indices_used)

  end subroutine test_scaling_performance_benchmark


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

  !> Test: Full outlier detection workflow.
  !> This test checks the complete workflow of outlier detection, including scaling, RDI calculation,
  !> outlier identification, and empirical p-values, for a small example with multiple families.
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
    real(real64)   :: p_values(n_genes)

    integer(int32) :: f, j, idx
    integer(int32) :: i_out, i_norm
    real(real64), parameter :: offsets(genes_per_fam) = [ &
        0.10_real64, 0.30_real64, 0.50_real64, 0.80_real64, 1.10_real64, 20.0_real64 ]
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
        is_outlier, loess_x, loess_y, loess_n, p_values, ierr, &
        60.0_real64)

    call assert_equal_int(ierr, 0, "ierr should be 0 (test_detect_outliers_basic)")
    call assert_true(any(is_outlier), "At least one outlier detected")

    ! --------------------------------------------------------------------------
    ! P-VALUE SANITY CHECKS
    ! --------------------------------------------------------------------------

    ! 1) All p-values must be in (0, 1] for this test (no negative RDIs expected)
    do idx = 1, n_genes
      call assert_true(p_values(idx) > 0.0_real64, "p-values must be > 0")
      call assert_true(p_values(idx) <= 1.0_real64, "p-values must be <= 1")
    end do

    ! 2) For each family: the constructed outlier (j=6) should have a smaller p-value
    !    than a typical in-family point (j=1).
    do f = 1, n_families
      i_norm = (f-1_int32)*genes_per_fam + 1_int32  ! j=1
      i_out  = (f-1_int32)*genes_per_fam + 6_int32  ! j=6 (big offset 20.0)

      call assert_true(p_values(i_out) <= p_values(i_norm), "Outlier must have <= p-value than normal")
    end do

    ! 3) If a gene is marked as outlier, it should have relatively small p-value.
    !    (We do not enforce a fixed alpha here because selection is percentile-based.)
    do idx = 1, n_genes
      if (is_outlier(idx)) then
        call assert_true(p_values(idx) < 0.50_real64, "Outliers should have p-values below 0.5 in this construction")
      end if
    end do

  end subroutine test_detect_outliers_basic

  !> Test: Outlier detection in a simple RDI vector using perm for ordering.
  !> Checks that identify_outliers marks the top values as outliers for the 75th percentile.
  subroutine test_identify_outliers_basic()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    ! Unsorted RDI (gene order)
    rdi = [0.3_real64, 0.1_real64, 0.4_real64, 0.2_real64]

    ! clamped copy (no negatives here)
    sorted_rdi = rdi

    ! perm such that sorted_rdi(perm) is ascending:
    ! values ascending: 0.1 (idx2), 0.2 (idx4), 0.3 (idx1), 0.4 (idx3)
    perm = [2_int32, 4_int32, 1_int32, 3_int32]

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 75.0_real64)

    ! percentile 75%: idx = ceil(4*0.75)=3
    ! threshold = sorted_rdi(perm(3)) = sorted_rdi(1) = 0.3
    call assert_equal_real(threshold, 0.3_real64, 1d-12, "Threshold at 75th percentile must be 0.3")

    ! Outliers: rdi >= 0.3 and >0 => genes 1 (0.3) and 3 (0.4)
    call assert_true(is_outlier(1), "Gene 1 (0.3) is outlier at 75th percentile")
    call assert_true(is_outlier(3), "Gene 3 (0.4) is outlier at 75th percentile")
    call assert_false(is_outlier(2), "Gene 2 (0.1) is not outlier")
    call assert_false(is_outlier(4), "Gene 4 (0.2) is not outlier")

    ! Validate p-values (upper-tail empirical with c=1), returned in gene order (same as rdi):
    ! distribution values = [0.1,0.2,0.3,0.4]
    ! d=0.3 => count_ge=2 (0.3,0.4) => p=(2+1)/5=0.6
    ! d=0.1 => count_ge=4 => p=(4+1)/5=1.0
    ! d=0.4 => count_ge=1 => p=(1+1)/5=0.4
    ! d=0.2 => count_ge=3 => p=(3+1)/5=0.8
    call assert_equal_real(p_values(1), 3.0_real64/denom, 1d-12, "p(gene1,d=0.3)=(2+c)/(n+c)=0.6")
    call assert_equal_real(p_values(2), 5.0_real64/denom, 1d-12, "p(gene2,d=0.1)=(4+c)/(n+c)=1.0")
    call assert_equal_real(p_values(3), 2.0_real64/denom, 1d-12, "p(gene3,d=0.4)=(1+c)/(n+c)=0.4")
    call assert_equal_real(p_values(4), 4.0_real64/denom, 1d-12, "p(gene4,d=0.2)=(3+c)/(n+c)=0.8")
  end subroutine test_identify_outliers_basic


  !> Test: All genes are outliers (percentile 0).
  !> At percentile=0, idx=1 => threshold = minimum. With rdi>0, all genes satisfy rdi>=min and rdi>0.
  subroutine test_all_outliers()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    rdi = [12.0_real64, 10.0_real64, 13.0_real64, 11.0_real64]
    sorted_rdi = rdi

    ! ascending by value: 10(idx2), 11(idx4), 12(idx1), 13(idx3)
    perm = [2_int32, 4_int32, 1_int32, 3_int32]

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 0.0_real64)

    call assert_equal_real(threshold, 10.0_real64, 1d-12, "Percentile 0 -> threshold must be min (10)")
    call assert_true(all(is_outlier), "All are outliers at 0 percentile (all rdi>0)")

    ! p-values sanity checks:
    ! distribution values = [10,11,12,13]
    ! gene2 has d=10 => count_ge=4 => p=(4+1)/5=1.0
    ! gene3 has d=13 => count_ge=1 => p=(1+1)/5=0.4
    call assert_equal_real(p_values(2), 5.0_real64/denom, 1d-12, "p(d=10)=1.0")
    call assert_equal_real(p_values(3), 2.0_real64/denom, 1d-12, "p(d=13)=0.4")
  end subroutine test_all_outliers


  !> Test: No genes are outliers (percentile 100).
  !> At percentile 100, threshold = maximum.
  !> With rule (rdi >= threshold), only the maximum (if positive) is outlier.
  subroutine test_no_outliers()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    sorted_rdi = rdi
    perm = [1_int32, 2_int32, 3_int32, 4_int32]

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 100.0_real64)

    call assert_equal_real(threshold, 0.4_real64, 1d-12, "Percentile 100 -> threshold must be maximum")
    call assert_false(any(is_outlier(1:3)), "No outliers below maximum at 100 percentile")
    call assert_true(is_outlier(4), "Maximum value must be outlier at 100 percentile")

    ! p-values: d=0.4 => count_ge=1 => p=(1+1)/5=0.4
    call assert_equal_real(p_values(4), 2.0_real64/denom, 1d-12, "p(max)=0.4 with c=1")
  end subroutine test_no_outliers


  !> Test: Percentile at boundaries (0 and 100).
  subroutine test_percentile_boundaries()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    real(real64) :: sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    sorted_rdi = rdi
    perm = [1_int32, 2_int32, 3_int32, 4_int32]

    ! Percentile 0 -> threshold = minimum (1.0)
    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 0.0_real64)
    call assert_equal_real(threshold, 1.0_real64, 1d-12, "p=0 -> threshold=min")
    call assert_true(all(is_outlier), "All positive values must be outliers at percentile 0")

    ! p-values spot checks (distribution [1,2,3,4]):
    ! d=1 => count_ge=4 => p=(4+1)/5=1.0
    ! d=4 => count_ge=1 => p=(1+1)/5=0.4
    call assert_equal_real(p_values(1), 5.0_real64/denom, 1d-12, "p(d=1)=1.0")
    call assert_equal_real(p_values(4), 2.0_real64/denom, 1d-12, "p(d=4)=0.4")

    ! Percentile 100 -> threshold = maximum (4.0)
    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 100.0_real64)
    call assert_equal_real(threshold, 4.0_real64, 1d-12, "p=100 -> threshold=max")
    call assert_false(any(is_outlier(1:3)), "Values below max not outliers at percentile 100")
    call assert_true(is_outlier(4), "Max value is outlier at percentile 100")
  end subroutine test_percentile_boundaries


  !> Test: Multiple genes with identical RDI at threshold.
  !> If multiple genes equal the threshold, all must be outliers.
  subroutine test_identical_rdi_at_threshold()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    sorted_rdi = rdi

    ! We need perm that makes sorted_rdi(perm) ascending.
    ! rdi already in nondecreasing order, perm identity is fine here.
    perm = [1_int32, 2_int32, 3_int32, 4_int32]

    ! percentile 50 -> idx = ceil(4*0.5)=2 -> threshold = second smallest = 2.0
    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 50.0_real64)

    call assert_equal_real(threshold, 2.0_real64, 1d-12, "Threshold must be 2.0 at 50 percentile")
    call assert_true(is_outlier(2), "First 2.0 must be outlier")
    call assert_true(is_outlier(3), "Second 2.0 must be outlier")
    call assert_false(is_outlier(1), "1.0 must not be outlier")

    ! p-values for d=2: count_ge=3 (2,2,3) => p=(3+1)/5=0.8
    call assert_equal_real(p_values(2), 4.0_real64/denom, 1d-12, "p(d=2)=0.8")
    call assert_equal_real(p_values(3), 4.0_real64/denom, 1d-12, "p(d=2)=0.8 (duplicate)")
  end subroutine test_identical_rdi_at_threshold


  !> Test: identify_outliers uses default percentile (95) if not provided.
  subroutine test_identify_outliers_default_percentile()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    real(real64), parameter :: c = 1.0_real64
    real(real64), parameter :: denom = real(n_genes, real64) + c  ! 5.0

    sorted_rdi = rdi
    perm = [1_int32, 2_int32, 3_int32, 4_int32]

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values)

    ! Default percentile 95: idx=ceil(4*0.95)=4 -> threshold=max=0.4
    call assert_equal_real(threshold, 0.4_real64, 1d-12, "Default percentile 95 -> threshold must be max")
    call assert_true(is_outlier(4), "Highest RDI is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")

    ! p(max)= (1+c)/(n+c)=2/5=0.4
    call assert_equal_real(p_values(4), 2.0_real64/denom, 1d-12, "p(max)=0.4")
  end subroutine test_identify_outliers_default_percentile


  !> Test: All RDI negative (should be ignored for outlier detection).
  !> This test checks that if all RDI values are negative, none are marked as outliers.
  subroutine test_all_negative_rdi()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=5
    real(real64) :: rdi(n_genes) = [-0.1_real64, -0.2_real64, -0.3_real64, -0.4_real64, -0.5_real64]
    real(real64) :: sorted_rdi(n_genes), p_values(n_genes)
    integer(int32) :: perm(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold

    ! Clamp negatives in the distribution to 0
    sorted_rdi = rdi
    where (sorted_rdi < 0.0_real64)
      sorted_rdi = 0.0_real64
    end where

    ! All values equal (0). Any perm is fine.
    perm = [1_int32, 2_int32, 3_int32, 4_int32, 5_int32]

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, 80.0_real64)

    call assert_equal_real(threshold, 0.0_real64, 1d-12, "All-negative -> clamped distribution -> threshold 0")
    call assert_true(.not. any(is_outlier), "All-negative RDI should not be outliers")

    ! Your compute_empirical_p_values rule: if d<0 => p=1 (per your earlier design),
    ! so all p-values must be 1 here.
    call assert_true(all(abs(p_values - 1.0_real64) < 1d-12), "All-negative rdi -> p_values must be 1.0 for all genes")
  end subroutine test_all_negative_rdi


  !> Exhaustive: detect_outliers error propagation (invalid indices).
  !> This test checks that detect_outliers propagates error_code 201 if gene_to_fam contains invalid indices.
  subroutine test_detect_outliers_invalid_indices()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,3,0]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: p_values(n_genes)
    integer(int32) :: ierr
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, p_values, ierr)
    call assert_equal_int(ierr, 201, 'detect_outliers propagates error_code 201 for invalid indices')
  end subroutine test_detect_outliers_invalid_indices

    !> Test: detect_outliers uses default percentile (95) if not provided.
  subroutine test_detect_outliers_default_percentile()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1_int32, 1_int32, 1_int32, 1_int32]

    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: ierr

    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    real(real64) :: p_values(n_genes)

    call detect_outliers( &
      n_genes, n_families, distances, gene_to_fam, &
      work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, p_values, ierr)

    call assert_equal_int(ierr, 0, "No error for default percentile")

    ! At 95th percentile, only the highest value should be outlier
    call assert_true(is_outlier(4), "Highest distance is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")

    ! p-value sanity (should be in (0,1] if computed)
    call assert_true(all(p_values > 0.0_real64), "p-values must be > 0")
    call assert_true(all(p_values <= 1.0_real64), "p-values must be <= 1")
  end subroutine test_detect_outliers_default_percentile


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
  !> This test checks that for a single gene in a single family, scaling and RDI are both 0.0,
  !> the gene is not an outlier, and the empirical p-value is 1.0.
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
    real(real64)   :: p_values(n_genes)

    ! 1) scaling (alloc wrapper)
    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, loess_n, ierr)
    call assert_equal_int(ierr, 0, "Scaling ierr should be 0 (test_single_gene_family)")

    ! 2) RDI
    perm = [(i, i=1,n_genes)]
    stack_left  = 0_int32
    stack_right = 0_int32

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! 3) Outlier detection (and p-values)
    perm = [(i, i=1,n_genes)]
    call detect_outliers( &
                        n_genes, n_families, distances, gene_to_fam, &
                        work_array, perm, stack_left, stack_right, &
                        is_outlier, loess_x, loess_y, loess_n, p_values, ierr, &
                        95.0_real64)

    call assert_equal_int(ierr, 0, "detect_outliers ierr should be 0 (test_single_gene_family)")

    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, "Single gene scaling (should be 0.0)")
    call assert_equal_real(rdi(1),    0.0_real64, 1e-12_real64, "Single gene RDI")
    call assert_false(is_outlier(1), "Single gene not outlier")

    ! Empirical p-value sanity for n=1: should be exactly 1.0 with c=1 and rdi=0
    call assert_equal_real(p_values(1), 1.0_real64, 1e-12_real64, "Single gene p-value should be 1.0")
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

  !> Test: Zero variance across families.
  !> Validates that the filter (is_close) handles identical log-SD values without excluding all data.
  subroutine test_outliers_zero_variance()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_families = 4_int32
    integer(int32), parameter :: genes_per_fam = 3_int32
    integer(int32), parameter :: n_genes = n_families * genes_per_fam

    real(real64)    :: distances(n_genes)
    integer(int32)  :: gene_to_fam(n_genes)
    real(real64)    :: dscale(n_families), loess_x(n_families), loess_y(n_families)
    integer(int32)  :: loess_n(n_families), ierr, i

    ! Setup: All genes have identical distances, resulting in identical family SDs
    distances = 1.0_real64
    do i = 1, n_genes
      gene_to_fam(i) = (i - 1) / genes_per_fam + 1
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, loess_n, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 (test_outliers_zero_variance)')
    ! Fallback should kick in, setting dscale to the median (0.0 in this specific case)
    call assert_true(all(dscale >= 0.0_real64), 'Dscale should be valid even with zero variance')
  end subroutine test_outliers_zero_variance

  !> Test: Orphan genes (Family Size = 1) with enough support for LOESS.
  subroutine test_outliers_orphan_genes()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    ! Usamos 8 familias para asegurar que n_eff >= 5
    integer(int32), parameter :: n_families = 8_int32
    integer(int32), parameter :: genes_per_fam = 3_int32
    integer(int32), parameter :: n_genes = (n_families-1)*genes_per_fam + 1

    real(real64)    :: distances(n_genes)
    integer(int32)  :: gene_to_fam(n_genes)
    real(real64)    :: dscale(n_families), lx(n_families), ly(n_families)
    integer(int32)  :: ln(n_families), ierr, i, j, idx

    ! Fill 7 families with valid data (3 genes each)
    idx = 0
    do i = 1, n_families - 1
      do j = 1, genes_per_fam
        idx = idx + 1
        gene_to_fam(idx) = i
        distances(idx) = real(i, real64) + (real(j, real64)*0.1)
      end do
    end do

    ! Family 8 is the Orphan (1 gene)
    idx = idx + 1
    gene_to_fam(idx) = 8
    distances(idx) = 10.0_real64

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      lx, ly, ln, ierr)

    call assert_equal_int(ierr, 0, 'Error code 0 (test_outliers_orphan_genes)')
    
    ! Check that valid families got a scaling factor
    call assert_true(dscale(1) > 0.0_real64, 'Valid family should have dscale > 0')
    
    ! Check that the orphan family (8) has zero scaling
    call assert_true(abs(dscale(8)) < 1e-12_real64, 'Orphan family should have dscale = 0')
  end subroutine

  !> Test: Extreme Outlier Detection.
  !> Verifies that a gene with a high RDI is correctly flagged as an outlier.
  subroutine test_outliers_extreme_detection()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), parameter :: n_genes = 12_int32
    integer(int32), parameter :: n_families = 2_int32

    integer(int32) :: gene_to_fam(12) = [ &
      1_int32,1_int32,1_int32,1_int32,1_int32,1_int32, &
      2_int32,2_int32,2_int32,2_int32,2_int32,2_int32 ]
    real(real64)    :: distances(n_genes)
    real(real64)    :: work_rdi(n_genes), lx(n_families), ly(n_families)
    integer(int32)  :: perm(n_genes), sl(n_genes), sr(n_genes), ln(n_families), ierr, i
    logical         :: is_outlier(n_genes)
    real(real64)    :: p_values(n_genes)

    do i = 1, n_genes - 1
        distances(i) = 0.8_real64 + (real(i, real64) * 0.05_real64)
    end do
    distances(n_genes) = 500.0_real64

    ! percentile=90: threshold at the 90th percentile; with n=12, this typically isolates the largest value
    perm = [(i, i=1,n_genes)]
    call detect_outliers( &
                         n_genes, n_families, distances, gene_to_fam, work_rdi, &
                         perm, sl, sr, is_outlier, lx, ly, ln, p_values, ierr, &
                         percentile=90.0_real64)

    call assert_equal_int(ierr, 0, "Error code 0 (test_outliers_extreme_detection)")
    call assert_true(is_outlier(12), "The extreme gene (12) should be an outlier")
    call assert_false(is_outlier(1), "The 1st gene should not be an outlier")

    ! p-value sanity: extreme should have smaller (or equal) p-value than a typical gene
    call assert_true(p_values(12) <= p_values(1), "Extreme outlier should have <= p-value than a normal gene")
    call assert_true(p_values(12) > 0.0_real64 .and. p_values(12) <= 1.0_real64, "p-values must be in (0,1]")
  end subroutine test_outliers_extreme_detection


  !> Test: NaN Handling.
  !> Validates that the presence of a NaN distance does not cause a crash and is handled safely.
  subroutine test_outliers_nan_handling()
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none

    integer(int32), parameter :: n_genes = 4_int32
    integer(int32), parameter :: n_families = 1_int32

    real(real64)    :: distances(n_genes)
    integer(int32)  :: gene_to_fam(n_genes) = [1_int32, 1_int32, 1_int32, 1_int32]
    real(real64)    :: work_rdi(n_genes), lx(n_families), ly(n_families)
    integer(int32)  :: perm(n_genes), sl(n_genes), sr(n_genes), ln(n_families), ierr, i
    logical         :: is_outlier(n_genes)
    real(real64)    :: p_values(n_genes)

    distances = 1.0_real64
    distances(1) = ieee_value(1.0_real64, ieee_quiet_nan)

    perm = [(i, i=1,n_genes)]
    call detect_outliers( &
                         n_genes, n_families, distances, gene_to_fam, work_rdi, &
                         perm, sl, sr, is_outlier, lx, ly, ln, p_values, ierr)

    call assert_equal_int(ierr, 0, "Should handle NaN without ierr (test_outliers_nan_handling)")

    ! Our design goal: NaNs should not be flagged as outliers.
    call assert_false(is_outlier(1), "NaN distance gene should not be marked as outlier")

    ! p-value sanity: should be in (0,1] (implementation may set NaN RDI to p=1 via guard logic)
    call assert_true(p_values(1) > 0.0_real64 .and. p_values(1) <= 1.0_real64, "NaN gene p-value must be in (0,1]")
  end subroutine test_outliers_nan_handling


end module mod_test_get_outliers


