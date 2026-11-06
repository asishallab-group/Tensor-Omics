!> The `group_centroid` cases, for all four entry points: `group_centroid_all` and
!| `group_centroid_orthologs`, and their `_expert` variants, which take the work array from the
!| caller. The centroid of each family, over all its genes or its orthologs only; unassigned
!| genes; empty families and their documented zero vector; gene order; cancellation and extreme magnitudes; the dimension and
!| range checks and the NaN/Inf check.
!|
!| Every value case runs through both variants of its mode (`check_all`, `check_orthologs`),
!| and every error case through all four (`check_error`). Every family but those of
!| `test_group_centroid_twenty_genes_per_family` has 1, 2 or 4 genes and exact sums, so its
!| centroid is exact even where an optimizer multiplies by the reciprocal of the count: those
!| cases need no tolerance.
module mod_test_group_centroid
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_gene_centroids, only: group_centroid_all, group_centroid_all_expert, &
                                  group_centroid_orthologs, group_centroid_orthologs_expert
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_group_centroid

    !> What `centroid_matrix` holds before a call, so that an output the procedure never writes
    !| shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64
    !> What the work array holds before an `_expert` call: every gene selected, so that a work
    !| array the procedure does not fully reset would add foreign genes to a centroid.
    logical(c_bool), parameter :: UNWRITTEN_MASK = .true._c_bool
    !> `gene_to_family` of a gene in no family (M_GENE_TO_FAM_SENTINEL in src/macros.h).
    integer(int32), parameter :: UNASSIGNED = 0_int32

contains

    !> Every case of `group_centroid_all`, `group_centroid_orthologs` and their `_expert` variants.
    function get_all_tests_group_centroid() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(12))
        all_tests(1) = test_case("test_group_centroid_all_values", test_group_centroid_all_values)
        all_tests(2) = test_case("test_group_centroid_orthologs_values", test_group_centroid_orthologs_values)
        all_tests(3) = test_case("test_group_centroid_unassigned_genes", test_group_centroid_unassigned_genes)
        all_tests(4) = test_case("test_group_centroid_empty_family", test_group_centroid_empty_family)
        all_tests(5) = test_case("test_group_centroid_gene_order_invariance", test_group_centroid_gene_order_invariance)
        all_tests(6) = test_case("test_group_centroid_twenty_genes_per_family", test_group_centroid_twenty_genes_per_family)
        all_tests(7) = test_case("test_group_centroid_single_gene", test_group_centroid_single_gene)
        all_tests(8) = test_case("test_group_centroid_cancellation", test_group_centroid_cancellation)
        all_tests(9) = test_case("test_group_centroid_extreme_magnitudes", test_group_centroid_extreme_magnitudes)
        all_tests(10) = test_case("test_group_centroid_dimensions", test_group_centroid_dimensions)
        all_tests(11) = test_case("test_group_centroid_gene_to_family_bounds", test_group_centroid_gene_to_family_bounds)
        all_tests(12) = test_case("test_group_centroid_rejects_nan_and_inf", test_group_centroid_rejects_nan_and_inf)
    end function get_all_tests_group_centroid

    ! --------------------------------------------------------------------------------------------
    ! Cases
    ! --------------------------------------------------------------------------------------------

    !> Seven genes, three families, interleaved: gene_to_family = [2, 1, 3, 1, 2, 1, 1], so family
    !| 1 has genes 2, 4, 6, 7; family 2 genes 1, 5; family 3 gene 3. Per family:
    !| axis 1: (1 + 2 + 3 + 6)/4 = 3, (10 + 30)/2 = 20, -7.5;
    !| axis 2: (0.5 + 0.25 + 1 + 0.25)/4 = 0.5, (-1 - 3)/2 = -2, 100.
    !| Families 1 and n_families (= 3) are the bounds of gene_to_family, both valid.
    subroutine test_group_centroid_all_values()
        integer(int32), parameter :: n_axes = 2, n_genes = 7, n_families = 3
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes)

        call seven_genes(vectors, gene_to_family)
        expected(:, 1) = [3.0_real64, 0.5_real64]
        expected(:, 2) = [20.0_real64, -2.0_real64]
        expected(:, 3) = [-7.5_real64, 100.0_real64]

        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, 0.0_real64, &
                       "test_group_centroid_all_values")
    end subroutine test_group_centroid_all_values

    !> The genes of `test_group_centroid_all_values`, with genes 2, 4 and 5 not orthologs. Each
    !| exclusion changes its family's centroid, so a filter that did nothing would show:
    !| family 1 keeps genes 6, 7: (3 + 6)/2 = 4.5, (1 + 0.25)/2 = 0.625;
    !| family 2 keeps gene 1: [10, -1]; family 3 keeps gene 3: [-7.5, 100].
    subroutine test_group_centroid_orthologs_values()
        integer(int32), parameter :: n_axes = 2, n_genes = 7, n_families = 3
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes)
        logical(c_bool) :: ortholog_set(n_genes)

        call seven_genes(vectors, gene_to_family)
        ortholog_set = [.true., .false., .true., .false., .false., .true., .true.]
        expected(:, 1) = [4.5_real64, 0.625_real64]
        expected(:, 2) = [10.0_real64, -1.0_real64]
        expected(:, 3) = [-7.5_real64, 100.0_real64]

        call check_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, expected, &
                             0.0_real64, "test_group_centroid_orthologs_values")
    end subroutine test_group_centroid_orthologs_values

    !> A gene in no family (gene_to_family = 0, the sentinel the wrapper accepts besides
    !| 1..n_families) is in no centroid, not even when it is an ortholog. Genes 2 and 4 are
    !| unassigned and hold 1000; family 1 is genes 1, 3: (2 + 4)/2 = 3, (-1 + 1)/2 = 0.
    subroutine test_group_centroid_unassigned_genes()
        integer(int32), parameter :: n_axes = 2, n_genes = 4, n_families = 1
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes)
        logical(c_bool) :: ortholog_set(n_genes)

        vectors(:, 1) = [2.0_real64, -1.0_real64]
        vectors(:, 2) = 1000.0_real64
        vectors(:, 3) = [4.0_real64, 1.0_real64]
        vectors(:, 4) = 1000.0_real64
        gene_to_family = [1, UNASSIGNED, 1, UNASSIGNED]
        ortholog_set = .true.
        expected(:, 1) = [3.0_real64, 0.0_real64]

        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, 0.0_real64, &
                       "test_group_centroid_unassigned_genes")
        call check_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, expected, &
                             0.0_real64, "test_group_centroid_unassigned_genes")
    end subroutine test_group_centroid_unassigned_genes

    !> Families without a gene to average: family 2 of three has no gene at all, and family 3,
    !| the last, has only gene 4, which is not an ortholog. Gene 3 is unassigned.
    !| All mode: family 1 = genes 1, 2: (4 + 8)/2 = 6; family 2 empty; family 3 = gene 4: 5.
    !| Orthologs mode: family 1 as before; families 2 and 3 empty.
    subroutine test_group_centroid_empty_family()
        integer(int32), parameter :: n_axes = 2, n_genes = 4, n_families = 3
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes)
        logical(c_bool) :: ortholog_set(n_genes)

        vectors(:, 1) = 4.0_real64
        vectors(:, 2) = 8.0_real64
        vectors(:, 3) = 1000.0_real64
        vectors(:, 4) = 5.0_real64
        gene_to_family = [1, 1, UNASSIGNED, 3]
        ortholog_set = [.true., .true., .true., .false.]

        ! A family with no gene to average gets the zero vector, with ierr = ERR_OK, as the doc of
        ! group_centroid_impl states: family 2 has no gene, and in the orthologs mode family 3 has
        ! none that is an ortholog.
        expected(:, 1) = 6.0_real64
        expected(:, 2) = 0.0_real64
        expected(:, 3) = 5.0_real64
        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, 0.0_real64, &
                       "test_group_centroid_empty_family")

        expected(:, 3) = 0.0_real64
        call check_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, expected, &
                             0.0_real64, "test_group_centroid_empty_family")
    end subroutine test_group_centroid_empty_family

    !> The genes of `test_group_centroid_orthologs_values` in another order, each carrying its
    !| vector, family and ortholog flag: [7, 3, 5, 1, 6, 2, 4]. A centroid is a mean over a set,
    !| so the expected values are those of the original order, in both modes.
    subroutine test_group_centroid_gene_order_invariance()
        integer(int32), parameter :: n_axes = 2, n_genes = 7, n_families = 3
        real(real64) :: vectors(n_axes, n_genes), shuffled_vectors(n_axes, n_genes)
        real(real64) :: expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes), shuffled_gene_to_family(n_genes), order(n_genes), i_gene
        logical(c_bool) :: ortholog_set(n_genes), shuffled_ortholog_set(n_genes)

        call seven_genes(vectors, gene_to_family)
        ortholog_set = [.true., .false., .true., .false., .false., .true., .true.]
        order = [7, 3, 5, 1, 6, 2, 4]
        do i_gene = 1, n_genes
            shuffled_vectors(:, i_gene) = vectors(:, order(i_gene))
            shuffled_gene_to_family(i_gene) = gene_to_family(order(i_gene))
            shuffled_ortholog_set(i_gene) = ortholog_set(order(i_gene))
        end do

        expected(:, 1) = [3.0_real64, 0.5_real64]
        expected(:, 2) = [20.0_real64, -2.0_real64]
        expected(:, 3) = [-7.5_real64, 100.0_real64]
        call check_all(shuffled_vectors, n_axes, n_genes, shuffled_gene_to_family, n_families, expected, &
                       0.0_real64, "test_group_centroid_gene_order_invariance")

        expected(:, 1) = [4.5_real64, 0.625_real64]
        expected(:, 2) = [10.0_real64, -1.0_real64]
        call check_orthologs(shuffled_vectors, n_axes, n_genes, shuffled_gene_to_family, n_families, &
                             shuffled_ortholog_set, expected, 0.0_real64, "test_group_centroid_gene_order_invariance")
    end subroutine test_group_centroid_gene_order_invariance

    !> A larger case, with a count that is not a power of two: 10 axes, 100 genes, 5 families,
    !| gene g in family mod(g - 1, 5) + 1 with the value g + 100*a on axis a. Family k holds the
    !| genes k + 5j, j = 0..19, so its sum on axis a is 20*(k + 100*a) + 5*(0 + ... + 19) =
    !| 20*(k + 100*a) + 950, and its centroid k + 100*a + 47.5.
    subroutine test_group_centroid_twenty_genes_per_family()
        integer(int32), parameter :: n_axes = 10, n_genes = 100, n_families = 5
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes), i_gene, i_axis, i_family

        do i_gene = 1, n_genes
            gene_to_family(i_gene) = mod(i_gene - 1, n_families) + 1
            do i_axis = 1, n_axes
                vectors(i_axis, i_gene) = real(i_gene + 100*i_axis, real64)
            end do
        end do
        do i_family = 1, n_families
            do i_axis = 1, n_axes
                expected(i_axis, i_family) = real(i_family + 100*i_axis, real64) + 47.5_real64
            end do
        end do

        ! The sums are exact integers and each division by 20 is exact as a division, but an
        ! optimizer may multiply by 1/20 instead, which is rounded and can cost a centroid its
        ! last bit: one ulp of the largest centroid, 5 + 1000 + 47.5 = 1052.5.
        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, spacing(1052.5_real64), &
                       "test_group_centroid_twenty_genes_per_family")
    end subroutine test_group_centroid_twenty_genes_per_family

    !> One gene in one family (n_genes = n_families = 1, the smallest valid sizes): its centroid
    !| is its vector, bit for bit, even for values no double holds exactly: x/1 is exact.
    !| n_axes = 1, the smallest valid size, is covered by a second call.
    subroutine test_group_centroid_single_gene()
        real(real64) :: vectors(3, 1), expected(3, 1), one_axis(1, 1)
        integer(int32) :: gene_to_family(1)
        logical(c_bool) :: ortholog_set(1)

        vectors(:, 1) = [12.3_real64, -4.5_real64, 6.7_real64]
        gene_to_family = [1]
        ortholog_set = .true.
        expected = vectors

        call check_all(vectors, 3, 1, gene_to_family, 1, expected, 0.0_real64, "test_group_centroid_single_gene")
        call check_orthologs(vectors, 3, 1, gene_to_family, 1, ortholog_set, expected, 0.0_real64, &
                             "test_group_centroid_single_gene")

        one_axis = -0.1_real64
        call check_all(one_axis, 1, 1, gene_to_family, 1, one_axis, 0.0_real64, &
                       "test_group_centroid_single_gene, one axis")
    end subroutine test_group_centroid_single_gene

    !> Values of very different magnitude that cancel: [2**40, -2**40, 0, 0] and
    !| [-2**-40, 2**-40, 5, -5] both average to 0. Every partial sum, in any order, fits in 53
    !| bits (5 - 2**-40 needs 43), so no rounding can leave a residue.
    subroutine test_group_centroid_cancellation()
        integer(int32), parameter :: n_axes = 2, n_genes = 4, n_families = 1
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families), big, small
        integer(int32) :: gene_to_family(n_genes)

        big = scale(1.0_real64, 40)
        small = scale(1.0_real64, -40)
        ! element by element: gfortran builds [big, -small] in a temporary
        vectors(1, 1) = big
        vectors(2, 1) = -small
        vectors(1, 2) = -big
        vectors(2, 2) = small
        vectors(:, 3) = [0.0_real64, 5.0_real64]
        vectors(:, 4) = [0.0_real64, -5.0_real64]
        gene_to_family = 1
        expected = 0.0_real64

        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, 0.0_real64, &
                       "test_group_centroid_cancellation")
    end subroutine test_group_centroid_cancellation

    !> A centroid is never larger than its largest value, so it must not overflow where the sum
    !| does. Per axis, over the two genes of the family: [huge, huge] -> huge;
    !| [2**1023, 1.5*2**1023] -> 1.25*2**1023; [tiny, 3*tiny] -> 2*tiny (normal numbers, no subnormal
    !| on the way); [huge, -huge] -> 0.
    subroutine test_group_centroid_extreme_magnitudes()
        integer(int32), parameter :: n_axes = 4, n_genes = 2, n_families = 1
        real(real64) :: vectors(n_axes, n_genes), expected(n_axes, n_families)
        integer(int32) :: gene_to_family(n_genes)
        logical(c_bool) :: ortholog_set(n_genes)

        vectors(:, 1) = [huge(1.0_real64), scale(1.0_real64, 1023), tiny(1.0_real64), huge(1.0_real64)]
        vectors(:, 2) = [huge(1.0_real64), scale(1.5_real64, 1023), 3.0_real64*tiny(1.0_real64), -huge(1.0_real64)]
        gene_to_family = 1
        ortholog_set = .true.
        expected(:, 1) = [huge(1.0_real64), scale(1.25_real64, 1023), 2.0_real64*tiny(1.0_real64), 0.0_real64]

        ! Regression: mean_vector_impl used to sum first and overflow: axes 1 and 2 came out as +Inf
        ! (huge + huge and 2**1023 + 1.5*2**1023 exceed huge, about 2**1024), an Inf from finite input.
        call check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, 0.0_real64, &
                       "test_group_centroid_extreme_magnitudes")
        call check_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, expected, &
                             0.0_real64, "test_group_centroid_extreme_magnitudes")
    end subroutine test_group_centroid_extreme_magnitudes

    !> n_axes (argument 2), n_genes (argument 3) and n_families (argument 5), each on its own:
    !| 0 is ERR_EMPTY_INPUT, -1 ERR_INVALID_INPUT. 1, the smallest valid size of each, is in
    !| `test_group_centroid_single_gene`.
    subroutine test_group_centroid_dimensions()
        real(real64) :: vectors(2, 2)
        integer(int32) :: gene_to_family(2)
        logical(c_bool) :: ortholog_set(2)

        vectors = 1.0_real64
        gene_to_family = [1, 1]
        ortholog_set = .true.

        call check_error(vectors, 0, 2, gene_to_family, 1, ortholog_set, ERR_EMPTY_INPUT, 2, &
                         "test_group_centroid_dimensions: n_axes = 0")
        call check_error(vectors, -1, 2, gene_to_family, 1, ortholog_set, ERR_INVALID_INPUT, 2, &
                         "test_group_centroid_dimensions: n_axes = -1")
        call check_error(vectors, 2, 0, gene_to_family, 1, ortholog_set, ERR_EMPTY_INPUT, 3, &
                         "test_group_centroid_dimensions: n_genes = 0")
        call check_error(vectors, 2, -1, gene_to_family, 1, ortholog_set, ERR_INVALID_INPUT, 3, &
                         "test_group_centroid_dimensions: n_genes = -1")
        ! gene_to_family is all unassigned here, so that no family index exceeds n_families < 1
        ! and blames argument 4 instead.
        gene_to_family = UNASSIGNED
        call check_error(vectors, 2, 2, gene_to_family, 0, ortholog_set, ERR_EMPTY_INPUT, 5, &
                         "test_group_centroid_dimensions: n_families = 0")
        call check_error(vectors, 2, 2, gene_to_family, -1, ortholog_set, ERR_INVALID_INPUT, 5, &
                         "test_group_centroid_dimensions: n_families = -1")
    end subroutine test_group_centroid_dimensions

    !> Every entry of gene_to_family (argument 4) is a family, 1..n_families, or the sentinel 0:
    !| -1 (just below the sentinel) and n_families + 1 are ERR_INVALID_INPUT. Each sits in the
    !| last entry, so the check must cover every gene. 0, 1 and n_families are valid in
    !| `test_group_centroid_all_values` and `test_group_centroid_unassigned_genes`.
    subroutine test_group_centroid_gene_to_family_bounds()
        integer(int32), parameter :: n_genes = 3, n_families = 2
        real(real64) :: vectors(1, n_genes)
        integer(int32) :: gene_to_family(n_genes)
        logical(c_bool) :: ortholog_set(n_genes)

        vectors = 1.0_real64
        ortholog_set = .true.

        gene_to_family = [1, 2, -1]
        call check_error(vectors, 1, n_genes, gene_to_family, n_families, ortholog_set, ERR_INVALID_INPUT, 4, &
                         "test_group_centroid_gene_to_family_bounds: -1")
        gene_to_family = [1, 2, n_families + 1]
        call check_error(vectors, 1, n_genes, gene_to_family, n_families, ortholog_set, ERR_INVALID_INPUT, 4, &
                         "test_group_centroid_gene_to_family_bounds: n_families + 1")
    end subroutine test_group_centroid_gene_to_family_bounds

    !> NaN, +Inf and -Inf in expression_vectors (argument 1) are ERR_NAN_INF. Each sits in the
    !| last element, of an unassigned gene: the whole matrix is checked, not only what is averaged.
    subroutine test_group_centroid_rejects_nan_and_inf()
        integer(int32), parameter :: n_axes = 2, n_genes = 2
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: vectors(n_axes, n_genes), bad(3)
        integer(int32) :: gene_to_family(n_genes), i_bad
        logical(c_bool) :: ortholog_set(n_genes)

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        vectors = 1.0_real64
        gene_to_family = [1, UNASSIGNED]
        ortholog_set = .true.
        do i_bad = 1, size(bad)
            vectors(n_axes, n_genes) = bad(i_bad)
            call check_error(vectors, n_axes, n_genes, gene_to_family, 1, ortholog_set, ERR_NAN_INF, 1, &
                             "test_group_centroid_rejects_nan_and_inf: "//trim(bad_names(i_bad)))
        end do
    end subroutine test_group_centroid_rejects_nan_and_inf

    ! --------------------------------------------------------------------------------------------
    ! Helpers
    ! --------------------------------------------------------------------------------------------

    !> The seven genes of the value cases: see `test_group_centroid_all_values`.
    subroutine seven_genes(vectors, gene_to_family)
        real(real64), intent(out) :: vectors(2, 7)
        integer(int32), intent(out) :: gene_to_family(7)

        vectors(:, 1) = [10.0_real64, -1.0_real64]
        vectors(:, 2) = [1.0_real64, 0.5_real64]
        vectors(:, 3) = [-7.5_real64, 100.0_real64]
        vectors(:, 4) = [2.0_real64, 0.25_real64]
        vectors(:, 5) = [30.0_real64, -3.0_real64]
        vectors(:, 6) = [3.0_real64, 1.0_real64]
        vectors(:, 7) = [6.0_real64, 0.25_real64]
        gene_to_family = [2, 1, 3, 1, 2, 1, 1]
    end subroutine seven_genes

    !> Runs `group_centroid_all` and `group_centroid_all_expert` on the same input, each into an
    !| output pre-filled with UNWRITTEN, and checks that both succeed with `expected`.
    subroutine check_all(vectors, n_axes, n_genes, gene_to_family, n_families, expected, tol, label)
        integer(int32), intent(in) :: n_axes, n_genes, n_families
        real(real64), intent(in) :: vectors(n_axes, n_genes), expected(n_axes, n_families), tol
        integer(int32), intent(in) :: gene_to_family(n_genes)
        character(*), intent(in) :: label
        real(real64) :: centroids(n_axes, n_families)
        logical(c_bool) :: work(n_genes)
        integer(int32) :: ierr

        centroids = UNWRITTEN
        call group_centroid_all(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, ierr)
        call assert_err(ierr, ERR_OK, label//": group_centroid_all, ierr")
        call assert_equal_array_real(centroids, expected, n_axes*n_families, tol, &
                                     label//": group_centroid_all, centroids", n_rows=n_axes)

        centroids = UNWRITTEN
        work = UNWRITTEN_MASK
        call group_centroid_all_expert(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, work, ierr)
        call assert_err(ierr, ERR_OK, label//": group_centroid_all_expert, ierr")
        call assert_equal_array_real(centroids, expected, n_axes*n_families, tol, &
                                     label//": group_centroid_all_expert, centroids", n_rows=n_axes)
    end subroutine check_all

    !> Runs `group_centroid_orthologs` and `group_centroid_orthologs_expert` on the same input,
    !| each into an output pre-filled with UNWRITTEN, and checks that both succeed with `expected`.
    subroutine check_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, expected, &
                               tol, label)
        integer(int32), intent(in) :: n_axes, n_genes, n_families
        real(real64), intent(in) :: vectors(n_axes, n_genes), expected(n_axes, n_families), tol
        integer(int32), intent(in) :: gene_to_family(n_genes)
        logical(c_bool), intent(in) :: ortholog_set(n_genes)
        character(*), intent(in) :: label
        real(real64) :: centroids(n_axes, n_families)
        logical(c_bool) :: work(n_genes)
        integer(int32) :: ierr

        centroids = UNWRITTEN
        call group_centroid_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, &
                                      ortholog_set, ierr)
        call assert_err(ierr, ERR_OK, label//": group_centroid_orthologs, ierr")
        call assert_equal_array_real(centroids, expected, n_axes*n_families, tol, &
                                     label//": group_centroid_orthologs, centroids", n_rows=n_axes)

        centroids = UNWRITTEN
        work = UNWRITTEN_MASK
        call group_centroid_orthologs_expert(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, &
                                             work, ortholog_set, ierr)
        call assert_err(ierr, ERR_OK, label//": group_centroid_orthologs_expert, ierr")
        call assert_equal_array_real(centroids, expected, n_axes*n_families, tol, &
                                     label//": group_centroid_orthologs_expert, centroids", n_rows=n_axes)
    end subroutine check_orthologs

    !> Runs all four entry points with dimensions that may be invalid, and checks that each
    !| rejects the input with `expected_code`, blaming argument `arg_pos` (the same position in
    !| all four). The arrays are as large as declared; the dimensions passed may claim less, or
    !| nothing, of them.
    subroutine check_error(vectors, n_axes, n_genes, gene_to_family, n_families, ortholog_set, &
                           expected_code, arg_pos, label)
        real(real64), contiguous, intent(in) :: vectors(:, :)
        integer(int32), intent(in) :: n_axes, n_genes, n_families, expected_code, arg_pos
        integer(int32), contiguous, intent(in) :: gene_to_family(:)
        logical(c_bool), contiguous, intent(in) :: ortholog_set(:)
        character(*), intent(in) :: label
        real(real64) :: centroids(size(vectors, 1), max(n_families, 1))
        logical(c_bool) :: work(size(vectors, 2))
        integer(int32) :: ierr

        call group_centroid_all(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, ierr)
        call assert_err(ierr, expected_code, label//": group_centroid_all", arg_pos)
        call group_centroid_all_expert(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, work, ierr)
        call assert_err(ierr, expected_code, label//": group_centroid_all_expert", arg_pos)
        call group_centroid_orthologs(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, &
                                      ortholog_set, ierr)
        call assert_err(ierr, expected_code, label//": group_centroid_orthologs", arg_pos)
        call group_centroid_orthologs_expert(vectors, n_axes, n_genes, gene_to_family, n_families, centroids, &
                                             work, ortholog_set, ierr)
        call assert_err(ierr, expected_code, label//": group_centroid_orthologs_expert", arg_pos)
    end subroutine check_error

end module mod_test_group_centroid
