! filepath: test/mod_test_euclidean_distance.f90
!> Unit test suite for euclidean distance routines.
module mod_test_euclidean_distance
  use asserts
  use tox_euclidean_distance
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(17)
    
    all_tests(1) = test_case("test_euclidean_distance_3d", test_euclidean_distance_3d)
    all_tests(2) = test_case("test_euclidean_distance_2d_origin", test_euclidean_distance_2d_origin)
    all_tests(3) = test_case("test_euclidean_distance_1d", test_euclidean_distance_1d)
    all_tests(4) = test_case("test_euclidean_distance_identical", test_euclidean_distance_identical)
    all_tests(5) = test_case("test_euclidean_distance_high_dim", test_euclidean_distance_high_dim)
    all_tests(6) = test_case("test_euclidean_distance_zero_vectors", test_euclidean_distance_zero_vectors)
    all_tests(7) = test_case("test_euclidean_distance_small_numbers", test_euclidean_distance_small_numbers)
    all_tests(8) = test_case("test_euclidean_distance_large_numbers", test_euclidean_distance_large_numbers)
    all_tests(9) = test_case("test_euclidean_distance_mixed_signs", test_euclidean_distance_mixed_signs)
    all_tests(10) = test_case("test_distance_to_centroid_basic", test_distance_to_centroid_basic)
    all_tests(11) = test_case("test_distance_to_centroid_single", test_distance_to_centroid_single)
    all_tests(12) = test_case("test_distance_to_centroid_high_dim", test_distance_to_centroid_high_dim)
    all_tests(13) = test_case("test_distance_to_centroid_invalid_families", test_distance_to_centroid_invalid_families)
    all_tests(14) = test_case("test_numerical_precision_epsilon", test_numerical_precision_epsilon)
    all_tests(15) = test_case("test_numerical_precision_close", test_numerical_precision_close)
    all_tests(16) = test_case("test_performance_large_scale", test_performance_large_scale)
    all_tests(17) = test_case("test_column_major_optimization", test_column_major_optimization)
  end function get_all_tests

  !> Run all euclidean distance tests.
  subroutine run_all_tests_euclidean_distance()
    type(test_case) :: all_tests(17)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All euclidean distance tests passed successfully."
  end subroutine run_all_tests_euclidean_distance

  !> Run specific euclidean distance tests by name.
  subroutine run_named_tests_euclidean_distance(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(17)
    integer(int32) :: i, j
    logical :: found
    
    all_tests = get_all_tests()
    
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
  end subroutine run_named_tests_euclidean_distance

  !> Test euclidean distance with standard 3D vectors.
  subroutine test_euclidean_distance_3d()
    real(real64) :: vec1(3), vec2(3), result, expected
    
    vec1 = [1.0_real64, 2.0_real64, 3.0_real64]
    vec2 = [4.0_real64, 5.0_real64, 6.0_real64]
    expected = sqrt(27.0_real64)
    
    call euclidean_distance(vec1, vec2, 3, result)
    call assert_equal_real(result, expected, 1e-12_real64, "Standard 3D vectors")
  end subroutine test_euclidean_distance_3d

  !> Test euclidean distance from 2D vector to origin.
  subroutine test_euclidean_distance_2d_origin()
    real(real64) :: vec1(2), vec2(2), result, expected
    
    vec1 = [3.0_real64, 4.0_real64]
    vec2 = [0.0_real64, 0.0_real64]
    expected = 5.0_real64
    
    call euclidean_distance(vec1, vec2, 2, result)
    call assert_equal_real(result, expected, 1e-12_real64, "2D vector to origin")
  end subroutine test_euclidean_distance_2d_origin

  !> Test euclidean distance with 1D vectors.
  subroutine test_euclidean_distance_1d()
    real(real64) :: vec1(1), vec2(1), result, expected
    
    vec1 = [5.0_real64]
    vec2 = [2.0_real64]
    expected = 3.0_real64
    
    call euclidean_distance(vec1, vec2, 1, result)
    call assert_equal_real(result, expected, 1e-12_real64, "1D vectors")
  end subroutine test_euclidean_distance_1d

  !> Test euclidean distance with identical vectors.
  subroutine test_euclidean_distance_identical()
    real(real64) :: vec1(3), vec2(3), result, expected
    
    vec1 = [1.0_real64, 2.0_real64, 3.0_real64]
    vec2 = [1.0_real64, 2.0_real64, 3.0_real64]
    expected = 0.0_real64
    
    call euclidean_distance(vec1, vec2, 3, result)
    call assert_equal_real(result, expected, 1e-15_real64, "Identical vectors")
  end subroutine test_euclidean_distance_identical

  !> Test euclidean distance with high-dimensional vectors.
  subroutine test_euclidean_distance_high_dim()
    real(real64) :: vec1(10), vec2(10), result, expected
    integer(int32) :: i
    
    ! Create test vectors with unit shift
    do i = 1, 10
      vec1(i) = real(i, real64)
      vec2(i) = real(i + 1, real64)
    end do
    expected = sqrt(10.0_real64)  ! sqrt(1^2 + 1^2 + ... + 1^2)
    
    call euclidean_distance(vec1, vec2, 10, result)
    call assert_equal_real(result, expected, 1e-12_real64, "10D unit-shift vectors")
  end subroutine test_euclidean_distance_high_dim

  !> Test euclidean distance with zero vectors.
  subroutine test_euclidean_distance_zero_vectors()
    real(real64) :: vec1(3), vec2(3), result, expected
    
    vec1 = [0.0_real64, 0.0_real64, 0.0_real64]
    vec2 = [0.0_real64, 0.0_real64, 0.0_real64]
    expected = 0.0_real64
    
    call euclidean_distance(vec1, vec2, 3, result)
    call assert_equal_real(result, expected, 1e-15_real64, "Zero vectors")
  end subroutine test_euclidean_distance_zero_vectors

  !> Test euclidean distance with very small numbers.
  subroutine test_euclidean_distance_small_numbers()
    real(real64) :: vec1(2), vec2(2), result, expected
    
    vec1 = [1e-15_real64, 1e-15_real64]
    vec2 = [2e-15_real64, 2e-15_real64]
    expected = sqrt(2.0_real64) * 1e-15_real64
    
    call euclidean_distance(vec1, vec2, 2, result)
    call assert_equal_real(result, expected, 1e-27_real64, "Very small numbers")
  end subroutine test_euclidean_distance_small_numbers

  !> Test euclidean distance with very large numbers.
  subroutine test_euclidean_distance_large_numbers()
    real(real64) :: vec1(2), vec2(2), result, expected
    
    vec1 = [1e15_real64, 1e15_real64]
    vec2 = [2e15_real64, 2e15_real64]
    expected = sqrt(2.0_real64) * 1e15_real64
    
    call euclidean_distance(vec1, vec2, 2, result)
    call assert_equal_real(result, expected, 1e3_real64, "Very large numbers")
  end subroutine test_euclidean_distance_large_numbers

  !> Test euclidean distance with mixed positive/negative values.
  subroutine test_euclidean_distance_mixed_signs()
    real(real64) :: vec1(2), vec2(2), result, expected
    
    vec1 = [-3.0_real64, 4.0_real64]
    vec2 = [3.0_real64, -4.0_real64]
    expected = 10.0_real64
    
    call euclidean_distance(vec1, vec2, 2, result)
    call assert_equal_real(result, expected, 1e-12_real64, "Mixed positive/negative")
  end subroutine test_euclidean_distance_mixed_signs

  !> Test basic distance_to_centroid functionality.
  subroutine test_distance_to_centroid_basic()
    integer(int32), parameter :: n_genes = 4, n_families = 2, d = 3
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    real(real64) :: expected_distances(n_genes)
    
    ! Gene data (column-major)
    genes(:, 1) = [1.0_real64, 0.0_real64, 0.0_real64]  ! Family 1
    genes(:, 2) = [0.0_real64, 1.0_real64, 0.0_real64]  ! Family 1
    genes(:, 3) = [3.0_real64, 0.0_real64, 0.0_real64]  ! Family 2
    genes(:, 4) = [0.0_real64, 3.0_real64, 0.0_real64]  ! Family 2
    
    ! Family assignments
    gene_to_fam = [1, 1, 2, 2]
    
    ! Centroids
    centroids(:, 1) = [0.5_real64, 0.5_real64, 0.0_real64]  ! Family 1 centroid
    centroids(:, 2) = [1.5_real64, 1.5_real64, 0.0_real64]  ! Family 2 centroid
    
    ! Expected distances
    expected_distances(1) = sqrt(0.5_real64**2 + 0.5_real64**2)  ! ~0.707
    expected_distances(2) = sqrt(0.5_real64**2 + 0.5_real64**2)  ! ~0.707
    expected_distances(3) = sqrt(1.5_real64**2 + 1.5_real64**2)  ! ~2.121
    expected_distances(4) = sqrt(1.5_real64**2 + 1.5_real64**2)  ! ~2.121
    
    ! Compute distances
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    ! Verify results
    call assert_equal_real(distances(1), expected_distances(1), 1e-12_real64, "Gene 1 distance")
    call assert_equal_real(distances(2), expected_distances(2), 1e-12_real64, "Gene 2 distance")
    call assert_equal_real(distances(3), expected_distances(3), 1e-12_real64, "Gene 3 distance")
    call assert_equal_real(distances(4), expected_distances(4), 1e-12_real64, "Gene 4 distance")
  end subroutine test_distance_to_centroid_basic

  !> Test distance_to_centroid with single gene/family.
  subroutine test_distance_to_centroid_single()
    integer(int32), parameter :: n_genes = 1, n_families = 1, d = 2
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    
    genes(:, 1) = [1.0_real64, 2.0_real64]
    centroids(:, 1) = [1.0_real64, 2.0_real64]
    gene_to_fam(1) = 1
    
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    call assert_equal_real(distances(1), 0.0_real64, 1e-15_real64, "Single gene/family")
  end subroutine test_distance_to_centroid_single

  !> Test distance_to_centroid with high-dimensional case.
  subroutine test_distance_to_centroid_high_dim()
    integer(int32), parameter :: n_genes = 3, n_families = 1, d = 100
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    integer(int32) :: i
    
    ! Initialize data
    do i = 1, d
      genes(i, 1) = real(i, real64)
      genes(i, 2) = real(i + 1, real64)
      genes(i, 3) = real(i + 2, real64)
      centroids(i, 1) = real(i + 1, real64)  ! Same as gene 2
    end do
    
    gene_to_fam = [1, 1, 1]
    
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    call assert_equal_real(distances(1), sqrt(real(d, real64)), 1e-10_real64, "High-dim gene 1")
    call assert_equal_real(distances(2), 0.0_real64, 1e-15_real64, "High-dim gene 2")
    call assert_equal_real(distances(3), sqrt(real(d, real64)), 1e-10_real64, "High-dim gene 3")
  end subroutine test_distance_to_centroid_high_dim

  !> Test distance_to_centroid with invalid family indices.
  subroutine test_distance_to_centroid_invalid_families()
    integer(int32), parameter :: n_genes = 3, n_families = 2, d = 2
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    
    genes(:, 1) = [1.0_real64, 2.0_real64]
    genes(:, 2) = [3.0_real64, 4.0_real64]
    genes(:, 3) = [5.0_real64, 6.0_real64]
    
    centroids(:, 1) = [0.0_real64, 0.0_real64]
    centroids(:, 2) = [1.0_real64, 1.0_real64]
    
    ! Invalid family indices
    gene_to_fam = [1, 3, 0]  ! family 3 doesn't exist, family 0 invalid
    
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    ! Check that valid gene has proper distance
    call assert_equal_real(distances(1), sqrt(5.0_real64), 1e-12_real64, "Valid gene distance")
    
    ! Check that invalid genes have error indicators
    call assert_equal_real(distances(2), -1.0_real64, 0.0_real64, "Invalid family index (3)")
    call assert_equal_real(distances(3), -1.0_real64, 0.0_real64, "Invalid family index (0)")
  end subroutine test_distance_to_centroid_invalid_families

  !> Test numerical precision near machine epsilon.
  subroutine test_numerical_precision_epsilon()
    real(real64), parameter :: eps = epsilon(1.0_real64)
    real(real64) :: vec1(2), vec2(2), result
    
    vec1 = [1.0_real64, 1.0_real64]
    vec2 = [1.0_real64 + eps, 1.0_real64 + eps]
    
    call euclidean_distance(vec1, vec2, 2, result)
    call assert_equal_real(result, sqrt(2.0_real64) * eps, eps * 10, "Machine epsilon precision")
  end subroutine test_numerical_precision_epsilon

  !> Test numerical precision with very close vectors.
  subroutine test_numerical_precision_close()
    real(real64) :: vec1(3), vec2(3), result
    
    vec1 = [1.234567890123456_real64, 2.345678901234567_real64, 3.456789012345678_real64]
    vec2 = [1.234567890123457_real64, 2.345678901234568_real64, 3.456789012345679_real64]
    
    call euclidean_distance(vec1, vec2, 3, result)
    call assert_equal_real(result, sqrt(3.0_real64) * 1e-15_real64, 1e-14_real64, "High precision vectors")
  end subroutine test_numerical_precision_close

  !> Test performance with large-scale data.
  subroutine test_performance_large_scale()
    integer(int32), parameter :: n_genes = 1000, n_families = 10, d = 50
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    integer(int32) :: i, j
    
    ! Initialize test data
    do i = 1, n_genes
      do j = 1, d
        genes(j, i) = sin(real(i * j, real64)) * cos(real(i + j, real64))
      end do
      gene_to_fam(i) = mod(i - 1, n_families) + 1
    end do
    
    do i = 1, n_families
      do j = 1, d
        centroids(j, i) = sin(real(i * j, real64)) * 0.5_real64
      end do
    end do
    
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    ! Verify some results are computed
    call assert_true(any(distances > 0.0_real64), "Performance test: distances computed")
  end subroutine test_performance_large_scale

  !> Test column-major memory optimization.
  subroutine test_column_major_optimization()
    integer(int32), parameter :: n_genes = 3, n_families = 2, d = 4
    real(real64) :: genes(d, n_genes), centroids(d, n_families)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: distances(n_genes)
    integer(int32) :: i
    
    ! Fill with known pattern to verify column-major access
    do i = 1, n_genes
      genes(:, i) = [real(i, real64), real(i*10, real64), &
                    real(i*100, real64), real(i*1000, real64)]
    end do
    
    centroids(:, 1) = [1.5_real64, 15.0_real64, 150.0_real64, 1500.0_real64]
    centroids(:, 2) = [2.5_real64, 25.0_real64, 250.0_real64, 2500.0_real64]
    
    gene_to_fam = [1, 1, 2]
    
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    ! Verify the pattern gives expected results
    ! Gene 1: [1,10,100,1000] vs Centroid 1: [1.5,15,150,1500]
    ! Distance should be sqrt(0.5^2 + 5^2 + 50^2 + 500^2) = sqrt(252525.25)
    call assert_equal_real(distances(1), sqrt(252525.25_real64), 1e-10_real64, "Column-major gene 1")
  end subroutine test_column_major_optimization

end module mod_test_euclidean_distance
