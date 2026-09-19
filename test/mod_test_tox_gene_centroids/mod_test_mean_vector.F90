!> The `mean_vector` cases: the element-wise mean of selected columns, the selection at its
!| bounds, the empty selection, cancellation, extreme magnitudes, the dimension and range checks
!| and the NaN/Inf check.
!|
!| Every selection but one has 1, 2 or 4 genes, and every sum is exact, so each mean is exact
!| even where an optimizer multiplies by the reciprocal of the count: those cases need no
!| tolerance. The one exception, `test_mean_vector_three_genes`, says why it has one.
module mod_test_mean_vector
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_gene_centroids, only: mean_vector
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_mean_vector

    !> What `centroid` holds before a call, so that an output the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64

contains

    !> Every case of `mean_vector`.
    function get_all_tests_mean_vector() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(12))
        all_tests(1) = test_case("test_mean_vector_values", test_mean_vector_values)
        all_tests(2) = test_case("test_mean_vector_selection_bounds", test_mean_vector_selection_bounds)
        all_tests(3) = test_case("test_mean_vector_three_genes", test_mean_vector_three_genes)
        all_tests(4) = test_case("test_mean_vector_no_genes", test_mean_vector_no_genes)
        all_tests(5) = test_case("test_mean_vector_repeated_index", test_mean_vector_repeated_index)
        all_tests(6) = test_case("test_mean_vector_cancellation", test_mean_vector_cancellation)
        all_tests(7) = test_case("test_mean_vector_extreme_magnitudes", test_mean_vector_extreme_magnitudes)
        all_tests(8) = test_case("test_mean_vector_dimensions", test_mean_vector_dimensions)
        all_tests(9) = test_case("test_mean_vector_n_selected_genes_bounds", test_mean_vector_n_selected_genes_bounds)
        all_tests(10) = test_case("test_mean_vector_gene_indices_bounds", test_mean_vector_gene_indices_bounds)
        all_tests(11) = test_case("test_mean_vector_rejects_nan_and_inf", test_mean_vector_rejects_nan_and_inf)
        all_tests(12) = test_case("test_mean_vector_single_axis_and_gene", test_mean_vector_single_axis_and_gene)
    end function get_all_tests_mean_vector

    !> Four of five genes, selected out of order: [5, 2, 1, 4]. Gene 3, left out, holds 1000 on
    !| every axis, so including it would show. Per axis, over genes 1, 2, 4, 5:
    !| (1 + 3 + 5 + 7)/4 = 4, (-2 + 6 - 10 + 2)/4 = -1, (8 + 0.5 + 2 + 1.5)/4 = 3.
    subroutine test_mean_vector_values()
        integer(int32), parameter :: n_axes = 3, n_genes = 5, n_selected = 4
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: gene_indices(n_selected), ierr

        vectors(:, 1) = [1.0_real64, -2.0_real64, 8.0_real64]
        vectors(:, 2) = [3.0_real64, 6.0_real64, 0.5_real64]
        vectors(:, 3) = 1000.0_real64
        vectors(:, 4) = [5.0_real64, -10.0_real64, 2.0_real64]
        vectors(:, 5) = [7.0_real64, 2.0_real64, 1.5_real64]
        gene_indices = [5, 2, 1, 4]
        expected = [4.0_real64, -1.0_real64, 3.0_real64]
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_selected, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_values: ierr")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, "test_mean_vector_values: centroid")
    end subroutine test_mean_vector_values

    !> The selection at its bounds: the first gene alone (index 1), the last alone (index
    !| n_genes), and every gene (n_selected_genes = n_genes, the largest valid count). A single
    !| gene's mean is that gene's vector, bit for bit, even for values like 12.3 that no double
    !| holds exactly: x/1 is exact.
    subroutine test_mean_vector_selection_bounds()
        integer(int32), parameter :: n_axes = 2, n_genes = 4
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: first(1), last(1), every(n_genes), ierr

        first = [1]
        last = [n_genes]
        every = [1, 2, 3, 4]

        vectors(:, 1) = [12.3_real64, -4.7_real64]
        vectors(:, 2) = [2.0_real64, 4.0_real64]
        vectors(:, 3) = [6.0_real64, -8.0_real64]
        vectors(:, 4) = [3.75_real64, 0.1_real64]

        expected = vectors(:, 1)
        centroid = UNWRITTEN
        call mean_vector(vectors, n_axes, n_genes, first, 1, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_selection_bounds: ierr, first gene")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, &
                                     "test_mean_vector_selection_bounds: centroid, first gene")

        expected = vectors(:, n_genes)
        centroid = UNWRITTEN
        call mean_vector(vectors, n_axes, n_genes, last, 1, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_selection_bounds: ierr, last gene")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, &
                                     "test_mean_vector_selection_bounds: centroid, last gene")

        ! Columns 1 and 4 replaced, so that the mean of all four is exact: over [0.25, -4.5],
        ! [2, 4], [6, -8], [3.75, 0.5], (0.25 + 2 + 6 + 3.75)/4 = 3 and (-4.5 + 4 - 8 + 0.5)/4 = -2.
        vectors(:, 1) = [0.25_real64, -4.5_real64]
        vectors(:, 4) = [3.75_real64, 0.5_real64]
        expected = [3.0_real64, -2.0_real64]
        centroid = UNWRITTEN
        call mean_vector(vectors, n_axes, n_genes, every, n_genes, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_selection_bounds: ierr, every gene")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, &
                                     "test_mean_vector_selection_bounds: centroid, every gene")
    end subroutine test_mean_vector_selection_bounds

    !> A count that is not a power of two: genes [1, 3, 5] per axis average to
    !| (1 + 3 + 5)/3 = 3, (2 - 4 + 8)/3 = 2 and (0.5 + 1 + 1.5)/3 = 1.
    subroutine test_mean_vector_three_genes()
        integer(int32), parameter :: n_axes = 3, n_genes = 3
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: gene_indices(n_genes), ierr

        vectors(:, 1) = [1.0_real64, 2.0_real64, 0.5_real64]
        vectors(:, 2) = [3.0_real64, -4.0_real64, 1.0_real64]
        vectors(:, 3) = [5.0_real64, 8.0_real64, 1.5_real64]
        gene_indices = [1, 2, 3]
        expected = [3.0_real64, 2.0_real64, 1.0_real64]
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_genes, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_three_genes: ierr")
        ! 9/3 is exact as a division, but an optimizer may multiply by 1/3 instead, which is
        ! rounded and can cost the result its last bit: one ulp of the largest value, 3.
        call assert_equal_array_real(centroid, expected, n_axes, spacing(3.0_real64), &
                                     "test_mean_vector_three_genes: centroid")
    end subroutine test_mean_vector_three_genes

    !> n_selected_genes = 0, the smallest valid count.
    subroutine test_mean_vector_no_genes()
        integer(int32), parameter :: n_axes = 3, n_genes = 2
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: gene_indices(0), ierr

        vectors = 5.0_real64
        ! QUESTION: the mean of no genes is undefined; today it is the zero vector, with ierr =
        ! ERR_OK. The published doc does not say so (only a comment in the implementation does),
        ! and a zero vector is indistinguishable from a gene expressed nowhere. Pinned as it is.
        expected = 0.0_real64
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, 0, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_no_genes: ierr")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, "test_mean_vector_no_genes: centroid")
    end subroutine test_mean_vector_no_genes

    !> An index selected more than once: [1, 1, 1, 2] over genes [0], [4], [100], [100]
    !| (n_selected_genes may not exceed n_genes, so the unselected genes 3 and 4 make room).
    subroutine test_mean_vector_repeated_index()
        integer(int32), parameter :: n_axes = 1, n_genes = 4, n_selected = 4
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: gene_indices(n_selected), ierr

        vectors(1, :) = [0.0_real64, 4.0_real64, 100.0_real64, 100.0_real64]
        gene_indices = [1, 1, 1, 2]
        ! QUESTION: the doc speaks of "the selected genes", a set, and caps their number at
        ! n_genes, but nothing rejects a repeated index. Today each occurrence counts:
        ! (0 + 0 + 0 + 4)/4 = 1, a weighted mean, not the mean of the two distinct genes,
        ! (0 + 4)/2 = 2. Pinned as it is.
        expected = [1.0_real64]
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_selected, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_repeated_index: ierr")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, &
                                     "test_mean_vector_repeated_index: centroid")
    end subroutine test_mean_vector_repeated_index

    !> Values of very different magnitude that cancel: [2**40, -2**40, 0, 0] and
    !| [-2**-40, 2**-40, 5, -5] both average to 0. Every partial sum, in any order, fits in 53
    !| bits (5 - 2**-40 needs 43), so no rounding can leave a residue.
    subroutine test_mean_vector_cancellation()
        integer(int32), parameter :: n_axes = 2, n_genes = 4
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes), big, small
        integer(int32) :: gene_indices(n_genes), ierr

        big = scale(1.0_real64, 40)
        small = scale(1.0_real64, -40)
        ! element by element: gfortran builds [big, -small] in a temporary
        vectors(1, 1) = big
        vectors(2, 1) = -small
        vectors(1, 2) = -big
        vectors(2, 2) = small
        vectors(:, 3) = [0.0_real64, 5.0_real64]
        vectors(:, 4) = [0.0_real64, -5.0_real64]
        gene_indices = [1, 2, 3, 4]
        expected = 0.0_real64
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_genes, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_cancellation: ierr")
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, "test_mean_vector_cancellation: centroid")
    end subroutine test_mean_vector_cancellation

    !> A mean is never larger than its largest value, so it must not overflow where the sum does.
    !| Per axis, over two genes: [huge, huge] -> huge; [2**1023, 1.5*2**1023] -> 1.25*2**1023;
    !| [tiny, 3*tiny] -> 2*tiny (normal numbers, no subnormal on the way); [huge, -huge] -> 0.
    subroutine test_mean_vector_extreme_magnitudes()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), expected(n_axes)
        integer(int32) :: gene_indices(n_genes), ierr

        vectors(:, 1) = [huge(1.0_real64), scale(1.0_real64, 1023), tiny(1.0_real64), huge(1.0_real64)]
        vectors(:, 2) = [huge(1.0_real64), scale(1.5_real64, 1023), 3.0_real64*tiny(1.0_real64), -huge(1.0_real64)]
        gene_indices = [1, 2]
        expected = [huge(1.0_real64), scale(1.25_real64, 1023), 2.0_real64*tiny(1.0_real64), 0.0_real64]
        centroid = UNWRITTEN

        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_genes, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_extreme_magnitudes: ierr")
        ! Regression: the sum used to come first and overflow: axes 1 and 2 came out as +Inf (huge +
        ! huge and 2**1023 + 1.5*2**1023 exceed huge, about 2**1024), an Inf from finite input.
        call assert_equal_array_real(centroid, expected, n_axes, 0.0_real64, &
                                     "test_mean_vector_extreme_magnitudes: centroid")
    end subroutine test_mean_vector_extreme_magnitudes

    !> n_axes (argument 2) and n_genes (argument 3), each on its own: 0 is ERR_EMPTY_INPUT, -1
    !| ERR_INVALID_INPUT. 1, the smallest valid size of both, is in
    !| `test_mean_vector_single_axis_and_gene`.
    subroutine test_mean_vector_dimensions()
        real(real64) :: vectors(2, 2), centroid(2)
        integer(int32) :: gene_indices(1), ierr

        vectors = 1.0_real64
        gene_indices = [1]

        call mean_vector(vectors, 0, 2, gene_indices, 1, centroid, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_mean_vector_dimensions: n_axes = 0", 2)
        call mean_vector(vectors, -1, 2, gene_indices, 1, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_dimensions: n_axes = -1", 2)
        ! n_selected_genes = 0 here, since with n_genes < 1 no count above 0 would be valid.
        call mean_vector(vectors, 2, 0, gene_indices, 0, centroid, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_mean_vector_dimensions: n_genes = 0", 3)
        call mean_vector(vectors, 2, -1, gene_indices, 0, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_dimensions: n_genes = -1", 3)
    end subroutine test_mean_vector_dimensions

    !> n_selected_genes (argument 5) is valid from 0 to n_genes: -1 and n_genes + 1 are
    !| ERR_INVALID_INPUT. 0 is in `test_mean_vector_no_genes`, n_genes in
    !| `test_mean_vector_selection_bounds`.
    subroutine test_mean_vector_n_selected_genes_bounds()
        integer(int32), parameter :: n_axes = 1, n_genes = 2
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes)
        integer(int32) :: gene_indices(n_genes + 1), ierr

        vectors = 1.0_real64
        gene_indices = [1, 2, 1]

        call mean_vector(vectors, n_axes, n_genes, gene_indices, -1, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_n_selected_genes_bounds: -1", 5)
        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_genes + 1, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_n_selected_genes_bounds: n_genes + 1", 5)
    end subroutine test_mean_vector_n_selected_genes_bounds

    !> Every entry of gene_indices (argument 4) is valid from 1 to n_genes: 0 and n_genes + 1
    !| are ERR_INVALID_INPUT. Each sits in the last entry, so the check must cover the whole
    !| selection. 1 and n_genes are in `test_mean_vector_selection_bounds`.
    subroutine test_mean_vector_gene_indices_bounds()
        integer(int32), parameter :: n_axes = 1, n_genes = 3, n_selected = 2
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes)
        integer(int32) :: gene_indices(n_selected), ierr

        vectors = 1.0_real64

        gene_indices = [1, 0]
        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_selected, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_gene_indices_bounds: index 0", 4)
        gene_indices = [1, n_genes + 1]
        call mean_vector(vectors, n_axes, n_genes, gene_indices, n_selected, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_mean_vector_gene_indices_bounds: index n_genes + 1", 4)
    end subroutine test_mean_vector_gene_indices_bounds

    !> NaN, +Inf and -Inf in expression_vectors (argument 1) are ERR_NAN_INF. Each sits in the
    !| last element, of a gene that is not even selected: the whole matrix is checked.
    subroutine test_mean_vector_rejects_nan_and_inf()
        integer(int32), parameter :: n_axes = 2, n_genes = 2
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: vectors(n_axes, n_genes), centroid(n_axes), bad(3)
        integer(int32) :: gene_indices(1), ierr, i_bad

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        vectors = 1.0_real64
        gene_indices = [1]
        do i_bad = 1, size(bad)
            vectors(n_axes, n_genes) = bad(i_bad)
            call mean_vector(vectors, n_axes, n_genes, gene_indices, 1, centroid, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_mean_vector_rejects_nan_and_inf: "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_mean_vector_rejects_nan_and_inf

    !> n_axes = n_genes = 1, the smallest valid sizes: the mean of the one gene is its value.
    subroutine test_mean_vector_single_axis_and_gene()
        real(real64) :: vectors(1, 1), centroid(1), expected(1)
        integer(int32) :: gene_indices(1), ierr

        vectors = -6.7_real64
        gene_indices = [1]
        expected = [-6.7_real64]
        centroid = UNWRITTEN

        call mean_vector(vectors, 1, 1, gene_indices, 1, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_mean_vector_single_axis_and_gene: ierr")
        call assert_equal_array_real(centroid, expected, 1, 0.0_real64, "test_mean_vector_single_axis_and_gene: centroid")
    end subroutine test_mean_vector_single_axis_and_gene

end module mod_test_mean_vector
