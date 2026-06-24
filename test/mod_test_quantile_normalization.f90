!> Unit test suite for quantile_normalization routine.
module mod_test_quantile_normalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use tox_normalization
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

    real(real64) :: TOL = 1d-12
    
contains

    !> Get array of all available tests.
    function get_all_tests_quantile_normalization() result(all_tests)
        type(test_case),allocatable :: all_tests(:)
        allocate(all_tests(9))

        all_tests(1) = test_case("test_error_zero_dimensions", test_error_zero_dimensions)
        all_tests(2) = test_case("test_error_negative_dimensions", test_error_negative_dimensions)
        all_tests(3) = test_case("test_trivial_1x1", test_trivial_1x1)
        all_tests(4) = test_case("test_single_tissue", test_single_tissue)
        all_tests(5) = test_case("test_single_gene", test_single_gene)
        all_tests(6) = test_case("test_simple_2x3", test_simple_2x3)
        all_tests(7) = test_case("test_ties", test_ties)
        all_tests(8) = test_case("test_already_normalized", test_already_normalized)
        all_tests(9) = test_case("test_random", test_random)
    end function get_all_tests_quantile_normalization

    ! ============================================================
    ! Small helper: stable, simple insertion sort for tests
    ! ============================================================
    subroutine sort_real_ascending(a, n)
        real(real64), intent(inout) :: a(n)
        integer(int32), intent(in)  :: n
        integer :: i, j
        real(real64) :: key
        do i = 2, n
            key = a(i)
            j = i - 1
            do while (j >= 1 .and. a(j) > key)
                a(j+1) = a(j)
                j = j - 1
            end do
            a(j+1) = key
        end do
    end subroutine


    ! ============================================================
    ! Test 1: Error handling — zero dimensions
    ! ============================================================
    subroutine test_error_zero_dimensions()
        real(real64) :: expr(1,1), norm(1,1), tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr = 1.0_real64

        call quantile_normalization(0, 1, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "n_genes=0 must trigger ERR_EMPTY_INPUT")

        call quantile_normalization(1, 0, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "n_tissues=0 must trigger ERR_EMPTY_INPUT")
    end subroutine


    ! ============================================================
    ! Test 2: Error handling — negative dimensions
    ! ============================================================
    subroutine test_error_negative_dimensions()
        real(real64) :: expr(1,1), norm(1,1), tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr = 1.0_real64

        call quantile_normalization(-3, 1, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "negative n_genes must trigger ERR_INVALID_INPUT")

        call quantile_normalization(1, -5, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "negative n_tissues must trigger ERR_INVALID_INPUT")
    end subroutine


    ! ============================================================
    ! Test 3: Trivial 1×1 case
    ! ============================================================
    subroutine test_trivial_1x1()
        real(real64) :: expr(1,1), norm(1,1), tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr(1,1) = 42.0_real64

        call quantile_normalization(1, 1, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "1x1: ierr must be OK")
        call assert_equal_real(norm(1,1), 42.0_real64, TOL, "1x1: value must remain unchanged")
        call assert_equal_real(means(1), 42.0_real64, TOL, "1x1: rank_means must equal value")
        call assert_permutation(perm, 1, "1x1: permutation must be [1]")
    end subroutine


    ! ============================================================
    ! Test 4: 1×N case (single tissue)
    ! ============================================================
    subroutine test_single_tissue()
        integer(int32), parameter :: n_genes=4, n_tissues=1
        real(real64) :: expr(1,n_genes), norm(1,n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr
        real(real64) :: sorted(n_genes)

        expr = reshape([3.0, 1.0, 4.0, 2.0], [1,n_genes])

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "1×N: ierr must be OK")

        sorted = expr(1,[2,4,1,3])
        call assert_equal_array_real(means, sorted, n_genes, TOL, "1×N: rank_means must equal sorted row")

        ! Each tissue must be permutation of rank_means
        call assert_equal_array_real(norm(1,[2,4,1,3]), means, n_genes, TOL, "1×N: result must match rank_means")
    end subroutine


    ! ============================================================
    ! Test 5: N×1 case (single gene)
    ! ============================================================
    subroutine test_single_gene()
        integer(int32), parameter :: n_genes=1, n_tissues=4
        real(real64) :: expr(n_tissues,1), norm(n_tissues,1)
        real(real64) :: tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr(:,1) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64]

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "N×1: ierr must be OK")

        ! All tissues get the same single rank mean
        call assert_equal_real(means(1), sum(expr(:,1))/n_tissues, TOL, "N×1: rank mean must be average")
        call assert_equal_array_real(norm(:,1), [means(1),means(1),means(1),means(1)], &
                                     n_tissues, TOL, "N×1: normalized column must be constant")
    end subroutine


    ! ============================================================
    ! Test 6: Simple 2×3 example with hand-computed expected values
    ! ============================================================
    subroutine test_simple_2x3()
        integer(int32), parameter :: n_genes=3, n_tissues=2
        real(real64) :: expr(n_tissues,n_genes), norm(n_tissues,n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr
        real(real64) :: s1(3), s2(3), expected_means(3)

        expr(1, :) = [3.0_real64,1.0_real64,2.0_real64]
        expr(2, :) = [6.0_real64,4.0_real64,5.0_real64]

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "2×3: ierr must be OK")

        expected_means = [2.5_real64, 3.5_real64, 4.5_real64]
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "2×3: rank_means correct")

        ! Check each tissue is permutation of rank_means
        s1 = norm(1,:)
        s2 = norm(2,:)
        call sort_real_ascending(s1, n_genes)
        call sort_real_ascending(s2, n_genes)

        call assert_equal_array_real(s1, expected_means, n_genes, TOL, "2×3: tissue 1 normalized correctly")
        call assert_equal_array_real(s2, expected_means, n_genes, TOL, "2×3: tissue 2 normalized correctly")
    end subroutine


    ! ============================================================
    ! Test 7: Ties
    ! ============================================================
    subroutine test_ties()
        integer(int32), parameter :: n_genes=4, n_tissues=2
        real(real64) :: expr(n_tissues,n_genes), norm(n_tissues,n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr
        real(real64) :: s1(4), s2(4)

        expr = reshape([1.0,2.0,2.0,3.0,  5.0,7.0,7.0,9.0], [n_tissues,n_genes])

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "ties: ierr must be OK")

        s1 = norm(1,:)
        s2 = norm(2,:)
        call sort_real_ascending(s1, n_genes)
        call sort_real_ascending(s2, n_genes)

        call assert_equal_array_real(s1, means, n_genes, TOL, "ties: tissue 1 sorted equals rank_means")
        call assert_equal_array_real(s2, means, n_genes, TOL, "ties: tissue 2 sorted equals rank_means")
    end subroutine


    ! ============================================================
    ! Test 8: Already normalized input
    ! ============================================================
    subroutine test_already_normalized()
        integer(int32), parameter :: n_genes=5, n_tissues=3
        real(real64) :: expr(n_tissues,n_genes), norm(n_tissues,n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr
        real(real64) :: sorted(n_genes)

        sorted = [1.0,2.0,3.0,4.0,5.0]
        expr(1,:) = sorted
        expr(2,:) = sorted
        expr(3,:) = sorted

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "already normalized: ierr must be OK")

        call assert_equal_array_real(means, sorted, n_genes, TOL, "already normalized: rank_means unchanged")
        call assert_equal_array_real(norm(1,:), sorted, n_genes, TOL, "already normalized: row 1 unchanged")
        call assert_equal_array_real(norm(2,:), sorted, n_genes, TOL, "already normalized: row 2 unchanged")
        call assert_equal_array_real(norm(3,:), sorted, n_genes, TOL, "already normalized: row 3 unchanged")
    end subroutine


    ! ============================================================
    ! Test 9: Random deterministic input
    ! ============================================================
    subroutine test_random()
        integer(int32), parameter :: n_genes=10, n_tissues=4
        real(real64) :: expr(n_tissues,n_genes), norm(n_tissues,n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr
        real(real64) :: row(n_genes)
        integer :: t,i

        ! Deterministic pseudo-random values
        do t = 1, n_tissues
            expr(t,:) = [( real(t*i,real64), i=1,n_genes )]
        end do

        call quantile_normalization(n_genes, n_tissues, expr, norm, means, tmp, perm, ierr)
        call assert_equal_int(ierr, ERR_OK, "random: ierr must be OK")

        call assert_sorted_real(means, n_genes, "random: rank_means must be sorted")

        do t = 1, n_tissues
            row = norm(t,:)
            call sort_real_ascending(row, n_genes)
            call assert_equal_array_real(row, means, n_genes, TOL, "random: each tissue matches rank_means")
        end do

        call assert_no_nan_real(means, n_genes, "random: no NaN in rank_means")
        call assert_no_inf_real(means, n_genes, "random: no Inf in rank_means")
    end subroutine

end module mod_test_quantile_normalization