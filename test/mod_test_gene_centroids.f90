! filepath: test/mod_test_gene_centroids.f90
!> Unit test suite for compute gene centroids routines.
module mod_test_gene_centroids
  use asserts
  use tox_gene_centroids
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT
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
    type(test_case) :: all_tests(10)
    all_tests(1) = test_case("test_basic_all_mode", test_basic_all_mode)
    all_tests(2) = test_case("test_basic_ortho_mode", test_basic_ortho_mode)
    all_tests(3) = test_case("test_empty_family", test_empty_family)
    all_tests(4) = test_case("test_no_matching_orthologs", test_no_matching_orthologs)
    all_tests(5) = test_case("test_single_gene_family", test_single_gene_family)
    all_tests(6) = test_case("test_extreme_values", test_extreme_values)
    all_tests(7) = test_case("test_higher_dimensions", test_higher_dimensions)
    all_tests(8) = test_case("test_gene_order_invariance", test_gene_order_invariance)
    all_tests(9) = test_case("test_invalid_input_arguments", test_invalid_input_arguments)
    all_tests(10) = test_case("test_invalid_family_mapping", test_invalid_family_mapping)
  end function get_all_tests

  !> Run all gene centroids tests.
  subroutine run_all_tests_gene_centroids()
    type(test_case) :: all_tests(10)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All gene centroids tests passed successfully."
  end subroutine run_all_tests_gene_centroids

  !> Run specific gene centroids tests by name.
  subroutine run_named_tests_gene_centroids(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(10)
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
  end subroutine run_named_tests_gene_centroids

  ! --------------------------------------------------------------------------
  ! Individual Test Cases
  ! --------------------------------------------------------------------------

  ! Test case 1: Basic functionality in 'all' mode.
  subroutine test_basic_all_mode()
    integer, parameter :: n_axes = 2, n_genes = 5, n_families = 2
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [n_axes, n_genes])
    gene_to_family = [1, 1, 2, 2, 1]
    ortholog_set = .false.
    expected = reshape([3.0, 3.0, 15.0, 15.0], [n_axes, n_families])

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, expected, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_basic_all_mode")
  end subroutine test_basic_all_mode

  ! Test case 2: Basic functionality in 'ortho' mode.
  subroutine test_basic_ortho_mode()
    integer, parameter :: n_axes = 2, n_genes = 5, n_families = 2
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [n_axes, n_genes])
    gene_to_family = [1, 1, 2, 2, 1]
    ortholog_set = [.true., .false., .true., .true., .true.]
    expected = reshape([3.0, 3.0, 15.0, 15.0], [n_axes, n_families])

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .false., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, expected, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_basic_ortho_mode")
  end subroutine test_basic_ortho_mode

  ! Test case 3: A family exists but has no genes assigned to it.
  subroutine test_empty_family()
    integer, parameter :: n_axes = 3, n_genes = 4, n_families = 2
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = 1.0
    gene_to_family = [1, 1, 1, 1] ! All genes in family 1, none in family 2
    ortholog_set = .true.
    expected = 0.0
    expected(:, 1) = 1.0 ! Centroid of family 1 is just the vector itself

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, expected, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_empty_family")
  end subroutine test_empty_family

  ! Test case 4: 'ortho' mode is selected, but a family has no orthologs.
  subroutine test_no_matching_orthologs()
    integer, parameter :: n_axes = 2, n_genes = 3, n_families = 1
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = reshape([10.0, 10.0, 20.0, 20.0, 30.0, 30.0], [n_axes, n_genes])
    gene_to_family = [1, 1, 1]
    ortholog_set = .false. ! No genes are orthologs
    expected = 0.0 ! Expect a zero vector

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .false., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, expected, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_no_matching_orthologs")
  end subroutine test_no_matching_orthologs

  ! Test case 5: A family contains only a single gene.
  subroutine test_single_gene_family()
    integer, parameter :: n_axes = 3, n_genes = 1, n_families = 1
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)

    vectors(:, 1) = [12.3, -4.5, 6.7]
    gene_to_family = [1]
    ortholog_set = .true.

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, vectors, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_single_gene_family")
  end subroutine test_single_gene_family

  ! Test case 6: Input vectors with extreme values.
  subroutine test_extreme_values()
    integer, parameter :: n_axes = 2, n_genes = 4, n_families = 1
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors(:, 1) = [1.0e12, -1.0e-12]
    vectors(:, 2) = [-1.0e12, 1.0e-12]
    vectors(:, 3) = [0.0, 5.0]
    vectors(:, 4) = [0.0, -5.0]
    gene_to_family = [1, 1, 1, 1]
    ortholog_set = .true.
    expected(:, 1) = [0.0, 0.0]

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids, expected, n_axes*n_families, 0.0_real64, 1e-9_real64, "test_extreme_values")
  end subroutine test_extreme_values

  ! Test case 7: Higher dimensional data.
  subroutine test_higher_dimensions()
    integer, parameter :: n_axes = 10, n_genes = 100, n_families = 5
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    integer :: i
    real(real64) :: expected(n_axes, n_families)

    do i = 1, n_genes
      vectors(:, i) = real(i, real64)
      gene_to_family(i) = mod(i - 1, n_families) + 1
    end do
    ortholog_set = .true.

    expected = 0.0
    expected(:, 1) = sum(vectors(:, 1:n_genes:n_families), dim=2)/real(count(gene_to_family == 1), real64)

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_allclose_array_real(centroids(:, 1), expected(:, 1), n_axes, 0.0_real64, 1e-9_real64, "test_higher_dimensions")
  end subroutine test_higher_dimensions

  ! Test case 8: Ensure the result is invariant to the order of genes.
  subroutine test_gene_order_invariance()
    integer, parameter :: n_axes = 2, n_genes = 5, n_families = 2
    real(real64) :: vectors1(n_axes, n_genes), centroids1(n_axes, n_families)
    integer(int32) :: gene_to_family1(n_genes), selected_indices1(n_genes), ierr
    logical :: ortholog_set1(n_genes)

    real(real64) :: vectors2(n_axes, n_genes), centroids2(n_axes, n_families)
    integer(int32) :: gene_to_family2(n_genes), selected_indices2(n_genes)
    logical :: ortholog_set2(n_genes)

    ! Setup 1: Original order
    vectors1 = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [n_axes, n_genes])
    gene_to_family1 = [1, 1, 2, 2, 1]
    ortholog_set1 = [.true., .false., .true., .true., .true.]

    ! Setup 2: Shuffled order
    vectors2 = reshape([5.0, 5.0, 10.0, 10.0, 1.0, 1.0, 3.0, 3.0, 20.0, 20.0], [n_axes, n_genes])
    gene_to_family2 = [1, 2, 1, 1, 2]
    ortholog_set2 = [.true., .true., .true., .false., .true.]

    call group_centroid(vectors1, n_axes, n_genes, gene_to_family1, n_families, &
                        centroids1, .false., ortholog_set1, selected_indices1, ierr)
    call group_centroid(vectors2, n_axes, n_genes, gene_to_family2, n_families, &
                        centroids2, .false., ortholog_set2, selected_indices2, ierr)

    call assert_allclose_array_real(centroids1, centroids2, n_axes*n_families, &
                                    0.0_real64, 1e-9_real64, "test_gene_order_invariance")
  end subroutine test_gene_order_invariance

  ! Test case 9: Test for invalid input arguments.
  subroutine test_invalid_input_arguments()
    integer, parameter :: n_axes = 0, n_genes = 5, n_families = 2, n_axes_invalid = 0, n_genes_invalid = 0, n_families_invalid = 0
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [n_axes, n_genes])
    gene_to_family = [1, 1, 2, 2, 1]
    ortholog_set = .false.
    expected = reshape([3.0, 3.0, 15.0, 15.0], [n_axes, n_families])

    call group_centroid(vectors, n_axes_invalid, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "Invalid 0 dimension should return ERR_EMPTY_INPUT")
    call group_centroid(vectors, n_axes, n_genes_invalid, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "Invalid 0 n_genes should return ERR_EMPTY_INPUT")
    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families_invalid, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "Invalid 0 families should return ERR_EMPTY_INPUT")
  end subroutine test_invalid_input_arguments

  ! Test case 10: Test for invalid family mapping.
  subroutine test_invalid_family_mapping()
    integer, parameter :: n_axes = 2, n_genes = 5, n_families = 2
    real(real64) :: vectors(n_axes, n_genes), centroids(n_axes, n_families)
    integer(int32) :: gene_to_family(n_genes), selected_indices(n_genes), ierr
    logical :: ortholog_set(n_genes)
    real(real64) :: expected(n_axes, n_families)

    vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [n_axes, n_genes])
    gene_to_family = [1, 1, 2, 3, 1] ! Invalid family mapping since family 3 does not exist
    ortholog_set = .false.
    expected = reshape([3.0, 3.0, 15.0, 15.0], [n_axes, n_families])

    call group_centroid(vectors, n_axes, n_genes, gene_to_family, n_families, &
                        centroids, .true., ortholog_set, selected_indices, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Invalid family mapping should return ERR_INVALID_INPUT")
  end subroutine test_invalid_family_mapping

end module mod_test_gene_centroids
