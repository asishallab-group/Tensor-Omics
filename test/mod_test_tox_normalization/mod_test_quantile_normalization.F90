!> The `quantile_normalization` cases: hand-derived rank means and the exact matrix each replicate
!| becomes, ties, the degenerate shapes, and the input checks.
module mod_test_quantile_normalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_normalization
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

    real(real64), parameter :: TOL = 1.0e-12_real64

contains

    !> Get array of all available tests.
    function get_all_tests_quantile_normalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(9))

        all_tests(1) = test_case("test_quantile_different_rank_orders", test_quantile_different_rank_orders)
        all_tests(2) = test_case("test_quantile_proportional_replicates", test_quantile_proportional_replicates)
        all_tests(3) = test_case("test_quantile_already_normalized", test_quantile_already_normalized)
        all_tests(4) = test_case("test_ties_share_the_mean_of_their_rank_means", &
                                 test_ties_share_the_mean_of_their_rank_means)
        all_tests(5) = test_case("test_quantile_single_replicate", test_quantile_single_replicate)
        all_tests(6) = test_case("test_quantile_single_gene", test_quantile_single_gene)
        all_tests(7) = test_case("test_quantile_trivial_1x1", test_quantile_trivial_1x1)
        all_tests(8) = test_case("test_quantile_rejects_nan_and_inf", test_quantile_rejects_nan_and_inf)
        all_tests(9) = test_case("test_quantile_dimensions", test_quantile_dimensions)
    end function get_all_tests_quantile_normalization

    !> The essence of quantile normalization: replicates that rank the genes differently. The
    !| sorted replicates [2, 4, 6] and [3, 6, 9] give rank means [2.5, 5, 7.5]. Replicate 1 ranks the
    !| genes 1, 2, 3 and gets them in that order; replicate 2 ranks gene 2 lowest and gene 1 highest.
    subroutine test_quantile_different_rank_orders()
        integer(int32), parameter :: n_genes = 3, n_replicates = 2
        real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
        real(real64) :: means(n_genes), expected_means(n_genes), tmp(n_genes)
        integer(int32) :: perm(n_genes), ierr

        expr(1, :) = [2.0_real64, 4.0_real64, 6.0_real64]
        expr(2, :) = [9.0_real64, 3.0_real64, 6.0_real64]
        expected_means = [2.5_real64, 5.0_real64, 7.5_real64]
        expected(1, :) = [2.5_real64, 5.0_real64, 7.5_real64]
        expected(2, :) = [7.5_real64, 2.5_real64, 5.0_real64]

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_different_rank_orders: ierr")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "test_quantile_different_rank_orders: rank means")
        call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                     "test_quantile_different_rank_orders: each gene must get the mean of its rank")
    end subroutine test_quantile_different_rank_orders

    !> Replicate t is t times [1, ..., 10]: every replicate ranks the genes alike, so gene i gets
    !| the mean of t*i over t = 1..4, which is 2.5*i, in every replicate.
    subroutine test_quantile_proportional_replicates()
        integer(int32), parameter :: n_genes = 10, n_replicates = 4
        real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
        real(real64) :: means(n_genes), expected_means(n_genes), tmp(n_genes)
        integer(int32) :: perm(n_genes), ierr, i_replicate, i_gene

        do i_gene = 1, n_genes
            expected_means(i_gene) = 2.5_real64*real(i_gene, real64)
            do i_replicate = 1, n_replicates
                expr(i_replicate, i_gene) = real(i_replicate*i_gene, real64)
                expected(i_replicate, i_gene) = expected_means(i_gene)
            end do
        end do

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_proportional_replicates: ierr")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "test_quantile_proportional_replicates: rank means")
        call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                     "test_quantile_proportional_replicates: every replicate must become 2.5*i")
    end subroutine test_quantile_proportional_replicates

    !> Replicates that already share one distribution come back unchanged.
    subroutine test_quantile_already_normalized()
        integer(int32), parameter :: n_genes = 5, n_replicates = 3
        real(real64), dimension(n_replicates, n_genes) :: expr, normalized
        real(real64) :: means(n_genes), expected_means(n_genes), tmp(n_genes)
        integer(int32) :: perm(n_genes), ierr, i_replicate

        expected_means = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        do i_replicate = 1, n_replicates
            expr(i_replicate, :) = expected_means
        end do

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_already_normalized: ierr")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "test_quantile_already_normalized: rank means")
        call assert_equal_array_real(normalized, expr, n_genes*n_replicates, TOL, &
                                     "test_quantile_already_normalized: must stay unchanged")
    end subroutine test_quantile_already_normalized

    !> Tied values within a replicate share the mean of the rank means their ranks span, so equal
    !| inputs stay equal. The rank means themselves do not change: ties only decide who gets them.
    subroutine test_ties_share_the_mean_of_their_rank_means()
        integer(int32), parameter :: n_genes = 4, n_replicates = 2
        real(real64) :: expr(n_replicates, n_genes), normalized(n_replicates, n_genes)
        real(real64) :: expected(n_replicates, n_genes), expected_means(n_genes)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr

        ! Replicate 1 ties genes 2 and 3 at ranks 2 and 3. The rank means of the sorted replicates
        ! [1, 3, 3, 5] and [2, 4, 6, 8] are [1.5, 3.5, 4.5, 6.5], so genes 2 and 3 of replicate 1
        ! both get (3.5 + 4.5)/2 = 4; replicate 2 has no ties.
        expr(1, :) = [1.0_real64, 3.0_real64, 3.0_real64, 5.0_real64]
        expr(2, :) = [2.0_real64, 4.0_real64, 6.0_real64, 8.0_real64]
        expected_means = [1.5_real64, 3.5_real64, 4.5_real64, 6.5_real64]
        expected(1, :) = [1.5_real64, 4.0_real64, 4.0_real64, 6.5_real64]
        expected(2, :) = expected_means

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "one tie: ierr must be OK")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "one tie: rank means")
        call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                     "one tie: genes 2 and 3 of replicate 1 must share (3.5 + 4.5)/2")

        ! A constant replicate is one tie over every rank. The rank means of [4, 4, 4, 4] and
        ! [1, 2, 3, 6] are [2.5, 3, 3.5, 5], so every gene of replicate 1 gets their mean, 3.5.
        expr(1, :) = 4.0_real64
        expr(2, :) = [1.0_real64, 2.0_real64, 3.0_real64, 6.0_real64]
        expected_means = [2.5_real64, 3.0_real64, 3.5_real64, 5.0_real64]
        expected(1, :) = 3.5_real64
        expected(2, :) = expected_means

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "constant replicate: ierr must be OK")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "constant replicate: rank means")
        call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                     "constant replicate: every gene must get the mean of all rank means")
    end subroutine test_ties_share_the_mean_of_their_rank_means

    !> A single replicate is its own reference distribution: it comes back unchanged, and the rank
    !| means are its sorted values.
    subroutine test_quantile_single_replicate()
        integer(int32), parameter :: n_genes = 4, n_replicates = 1
        real(real64), dimension(n_replicates, n_genes) :: expr, normalized
        real(real64) :: means(n_genes), expected_means(n_genes), tmp(n_genes)
        integer(int32) :: perm(n_genes), ierr

        expr(1, :) = [3.0_real64, 1.0_real64, 4.0_real64, 2.0_real64]
        expected_means = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_single_replicate: ierr")
        call assert_equal_array_real(means, expected_means, n_genes, TOL, "test_quantile_single_replicate: rank means")
        call assert_equal_array_real(normalized, expr, n_genes, TOL, "test_quantile_single_replicate: must stay unchanged")
    end subroutine test_quantile_single_replicate

    !> A single gene has one rank: every replicate gets the mean of its values, here 25.
    subroutine test_quantile_single_gene()
        integer(int32), parameter :: n_genes = 1, n_replicates = 4
        real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
        real(real64) :: means(n_genes), tmp(n_genes)
        integer(int32) :: perm(n_genes), ierr

        expr(:, 1) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64]
        expected = 25.0_real64

        call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_single_gene: ierr")
        call assert_equal_real(means(1), 25.0_real64, TOL, "test_quantile_single_gene: rank mean")
        call assert_equal_array_real(normalized, expected, n_replicates, TOL, "test_quantile_single_gene: every replicate gets 25")
    end subroutine test_quantile_single_gene

    !> A 1 x 1 matrix is unchanged, and so is its one rank mean.
    subroutine test_quantile_trivial_1x1()
        real(real64) :: expr(1, 1), normalized(1, 1), tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr(1, 1) = 42.0_real64

        call quantile_normalization_expert(1, 1, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_quantile_trivial_1x1: ierr")
        call assert_equal_real(normalized(1, 1), 42.0_real64, 0.0_real64, "test_quantile_trivial_1x1: value unchanged")
        call assert_equal_real(means(1), 42.0_real64, 0.0_real64, "test_quantile_trivial_1x1: rank mean")
    end subroutine test_quantile_trivial_1x1

    !> NaN and Inf are rejected up front rather than carried into the result: TOX writes no NaN.
    !| Unchecked, one NaN would spread through the rank means into every replicate.
    subroutine test_quantile_rejects_nan_and_inf()
        integer(int32), parameter :: n_genes = 3, n_replicates = 2
        real(real64) :: expr(n_replicates, n_genes), normalized(n_replicates, n_genes), bad(2)
        real(real64) :: tmp(n_genes), means(n_genes)
        integer(int32) :: perm(n_genes), ierr, i_bad

        expr(1, :) = [1.0_real64, 2.0_real64, 3.0_real64]
        expr(2, :) = [4.0_real64, 5.0_real64, 6.0_real64]
        bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

        do i_bad = 1, size(bad)
            expr(1, 1) = bad(i_bad)
            call quantile_normalization_expert(n_genes, n_replicates, expr, normalized, means, tmp, perm, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                                  "test_quantile_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
        end do
    end subroutine test_quantile_rejects_nan_and_inf

    !> Each dimension on its own: zero is ERR_EMPTY_INPUT, negative ERR_INVALID_INPUT.
    subroutine test_quantile_dimensions()
        real(real64) :: expr(1, 1), normalized(1, 1), tmp(1), means(1)
        integer(int32) :: perm(1), ierr

        expr = 1.0_real64
        call quantile_normalization_expert(0, 1, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_quantile_dimensions: n_genes = 0")
        call quantile_normalization_expert(1, 0, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_quantile_dimensions: n_replicates = 0")
        call quantile_normalization_expert(-3, 1, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_quantile_dimensions: n_genes = -3")
        call quantile_normalization_expert(1, -5, expr, normalized, means, tmp, perm, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_quantile_dimensions: n_replicates = -5")
    end subroutine test_quantile_dimensions

end module mod_test_quantile_normalization
