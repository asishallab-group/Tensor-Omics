!> The `compute_shift_vector_field` cases: for every gene, its family's centroid and the shift
!| from that centroid to the gene, `expression - centroid`; genes sharing a family, one family
!| per gene, a gene at its centroid, unassigned genes (the sentinel 0, which has no centroid and
!| gives two zero vectors), more families than genes, many genes, extreme magnitudes, the
!| smallest sizes, the family-id range, the dimension checks and the NaN/Inf checks.
!|
!| Every expected value is written out by hand. The inputs are small integers, exact powers of
!| two or `huge` and its exact halves, so every difference is exact (or overflows by design, as
!| one case documents) and every comparison has tolerance 0.
module mod_test_compute_shift_vector_field
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_shift_vectors, only: compute_shift_vector_field
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_compute_shift_vector_field

    !> What `shift_vectors` holds before a call, so that an element the procedure never writes
    !| shows. No expected value in this suite equals it.
    real(real64), parameter :: UNWRITTEN = -777.0_real64
    !> `gene_to_fam` of a gene in no family (M_GENE_TO_FAM_SENTINEL in src/macros.h).
    integer(int32), parameter :: UNASSIGNED = 0_int32

    !> How many non-finite values `non_finite_value` hands out.
    integer(int32), parameter :: N_NON_FINITE = 3
    !> The names of the non-finite values, for the assertion messages.
    character(*), parameter :: NON_FINITE_NAMES(N_NON_FINITE) = ["NaN ", "+Inf", "-Inf"]

contains

    !> Every case of `compute_shift_vector_field`.
    function get_all_tests_compute_shift_vector_field() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(13))
        all_tests(1) = test_case("test_compute_shift_vector_field_values", test_compute_shift_vector_field_values)
        all_tests(2) = test_case("test_compute_shift_vector_field_shared_family", &
                                 test_compute_shift_vector_field_shared_family)
        all_tests(3) = test_case("test_compute_shift_vector_field_one_family_per_gene", &
                                 test_compute_shift_vector_field_one_family_per_gene)
        all_tests(4) = test_case("test_compute_shift_vector_field_gene_at_centroid", &
                                 test_compute_shift_vector_field_gene_at_centroid)
        all_tests(5) = test_case("test_compute_shift_vector_field_unassigned_genes", &
                                 test_compute_shift_vector_field_unassigned_genes)
        all_tests(6) = test_case("test_compute_shift_vector_field_all_unassigned", &
                                 test_compute_shift_vector_field_all_unassigned)
        all_tests(7) = test_case("test_compute_shift_vector_field_unused_families", &
                                 test_compute_shift_vector_field_unused_families)
        all_tests(8) = test_case("test_compute_shift_vector_field_many_genes", test_compute_shift_vector_field_many_genes)
        all_tests(9) = test_case("test_compute_shift_vector_field_extreme_magnitudes", &
                                 test_compute_shift_vector_field_extreme_magnitudes)
        all_tests(10) = test_case("test_compute_shift_vector_field_smallest_sizes", &
                                  test_compute_shift_vector_field_smallest_sizes)
        all_tests(11) = test_case("test_compute_shift_vector_field_family_id_range", &
                                  test_compute_shift_vector_field_family_id_range)
        all_tests(12) = test_case("test_compute_shift_vector_field_dimensions", test_compute_shift_vector_field_dimensions)
        all_tests(13) = test_case("test_compute_shift_vector_field_rejects_nan_and_inf", &
                                  test_compute_shift_vector_field_rejects_nan_and_inf)
    end function get_all_tests_compute_shift_vector_field

    !> Calls the procedure and checks `ierr` OK and every element of `shift_vectors` exactly.
    !| `shift_vectors` starts at UNWRITTEN and has one gene more than the call is told, so an
    !| element the procedure skips and a write past the last gene both show.
    subroutine check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                   expected, msg)
        integer(int32), intent(in) :: n_tissues, n_genes, n_families
        real(real64), intent(in) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
        real(real64), intent(in) :: expected(n_tissues, 2, n_genes)
        character(*), intent(in) :: msg
        real(real64) :: shift_vectors(n_tissues, 2, n_genes + 1), unwritten_gene(n_tissues, 2)
        integer(int32) :: ierr

        shift_vectors = UNWRITTEN
        unwritten_gene = UNWRITTEN
        call compute_shift_vector_field(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                        shift_vectors, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        ! flattened (tissue, vector, gene): row = tissue, column = 2*(gene - 1) + vector
        call assert_equal_array_real(shift_vectors, expected, n_tissues*2*n_genes, 0.0_real64, msg//": shift_vectors", &
                                     n_rows=n_tissues)
        call assert_equal_array_real(shift_vectors(:, :, n_genes + 1), unwritten_gene, n_tissues*2, 0.0_real64, &
                                     msg//": no write past the last gene")
    end subroutine check_shift_vectors

    !> The `i_value`-th non-finite value: NaN, +Inf, -Inf.
    function non_finite_value(i_value) result(value)
        integer(int32), intent(in) :: i_value
        real(real64) :: value

        select case (i_value)
        case (1)
            value = ieee_value(1.0_real64, ieee_quiet_nan)
        case (2)
            value = ieee_value(1.0_real64, ieee_positive_inf)
        case default
            value = ieee_value(1.0_real64, ieee_negative_inf)
        end select
    end function non_finite_value

    !> Kept from the old suite (`test_correct_family_mapping`): 3 tissues, 5 genes, 3 families,
    !| the genes in no order of their families, every family used, the first and the last id.
    !| Genes 1..5 are [1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]; families
    !| 1..3 are [5, 4, 3], [2, 1, 0], [-1, -2, -3]; gene_to_fam = [2, 3, 1, 3, 1]:
    !| - gene 1, family 2: shift [1-2, 2-1, 3-0] = [-1, 1, 3];
    !| - gene 2, family 3: [4+1, 5+2, 6+3] = [5, 7, 9];
    !| - gene 3, family 1: [7-5, 8-4, 9-3] = [2, 4, 6];
    !| - gene 4, family 3: [10+1, 11+2, 12+3] = [11, 13, 15];
    !| - gene 5, family 1: [13-5, 14-4, 15-3] = [8, 10, 12].
    subroutine test_compute_shift_vector_field_values()
        integer(int32), parameter :: n_tissues = 3, n_genes = 5, n_families = 3
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        expression_vectors(:, 2) = [4.0_real64, 5.0_real64, 6.0_real64]
        expression_vectors(:, 3) = [7.0_real64, 8.0_real64, 9.0_real64]
        expression_vectors(:, 4) = [10.0_real64, 11.0_real64, 12.0_real64]
        expression_vectors(:, 5) = [13.0_real64, 14.0_real64, 15.0_real64]
        family_centroids(:, 1) = [5.0_real64, 4.0_real64, 3.0_real64]
        family_centroids(:, 2) = [2.0_real64, 1.0_real64, 0.0_real64]
        family_centroids(:, 3) = [-1.0_real64, -2.0_real64, -3.0_real64]
        gene_to_fam = [2, 3, 1, 3, 1]

        expected(:, 1, 1) = [2.0_real64, 1.0_real64, 0.0_real64]
        expected(:, 2, 1) = [-1.0_real64, 1.0_real64, 3.0_real64]
        expected(:, 1, 2) = [-1.0_real64, -2.0_real64, -3.0_real64]
        expected(:, 2, 2) = [5.0_real64, 7.0_real64, 9.0_real64]
        expected(:, 1, 3) = [5.0_real64, 4.0_real64, 3.0_real64]
        expected(:, 2, 3) = [2.0_real64, 4.0_real64, 6.0_real64]
        expected(:, 1, 4) = [-1.0_real64, -2.0_real64, -3.0_real64]
        expected(:, 2, 4) = [11.0_real64, 13.0_real64, 15.0_real64]
        expected(:, 1, 5) = [5.0_real64, 4.0_real64, 3.0_real64]
        expected(:, 2, 5) = [8.0_real64, 10.0_real64, 12.0_real64]
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_values")
    end subroutine test_compute_shift_vector_field_values

    !> Kept from the old suite (`test_multiple_genes_per_family`), which compared 4 of the 16
    !| elements and expected wrong values in 5 of the 12 it never compared (3 vectors). 2 tissues, genes
    !| [1, 2], [3, 4], [5, 6], [7, 8] in families [1, 2, 1, 2], families [10, 20] and [30, 40]:
    !| - gene 1: centroid [10, 20], shift [-9, -18];
    !| - gene 2: centroid [30, 40], shift [-27, -36];
    !| - gene 3: centroid [10, 20], shift [-5, -14];
    !| - gene 4: centroid [30, 40], shift [-23, -32].
    subroutine test_compute_shift_vector_field_shared_family()
        integer(int32), parameter :: n_tissues = 2, n_genes = 4, n_families = 2
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64]
        expression_vectors(:, 2) = [3.0_real64, 4.0_real64]
        expression_vectors(:, 3) = [5.0_real64, 6.0_real64]
        expression_vectors(:, 4) = [7.0_real64, 8.0_real64]
        family_centroids(:, 1) = [10.0_real64, 20.0_real64]
        family_centroids(:, 2) = [30.0_real64, 40.0_real64]
        gene_to_fam = [1, 2, 1, 2]

        expected(:, 1, 1) = [10.0_real64, 20.0_real64]
        expected(:, 2, 1) = [-9.0_real64, -18.0_real64]
        expected(:, 1, 2) = [30.0_real64, 40.0_real64]
        expected(:, 2, 2) = [-27.0_real64, -36.0_real64]
        expected(:, 1, 3) = [10.0_real64, 20.0_real64]
        expected(:, 2, 3) = [-5.0_real64, -14.0_real64]
        expected(:, 1, 4) = [30.0_real64, 40.0_real64]
        expected(:, 2, 4) = [-23.0_real64, -32.0_real64]
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_shared_family")
    end subroutine test_compute_shift_vector_field_shared_family

    !> Kept from the old suite (`test_single_gene_per_family`), which compared 4 of the 16
    !| elements. Every gene its own family, in reverse order, so that gene i reading family i
    !| would show: genes [1, 2], [3, 4], [5, 6], [7, 8] in families [4, 3, 2, 1], families
    !| [10, 20], [30, 40], [50, 60], [70, 80]:
    !| - gene 1, family 4: centroid [70, 80], shift [-69, -78];
    !| - gene 2, family 3: centroid [50, 60], shift [-47, -56];
    !| - gene 3, family 2: centroid [30, 40], shift [-25, -34];
    !| - gene 4, family 1: centroid [10, 20], shift [-3, -12].
    subroutine test_compute_shift_vector_field_one_family_per_gene()
        integer(int32), parameter :: n_tissues = 2, n_genes = 4, n_families = 4
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64]
        expression_vectors(:, 2) = [3.0_real64, 4.0_real64]
        expression_vectors(:, 3) = [5.0_real64, 6.0_real64]
        expression_vectors(:, 4) = [7.0_real64, 8.0_real64]
        family_centroids(:, 1) = [10.0_real64, 20.0_real64]
        family_centroids(:, 2) = [30.0_real64, 40.0_real64]
        family_centroids(:, 3) = [50.0_real64, 60.0_real64]
        family_centroids(:, 4) = [70.0_real64, 80.0_real64]
        gene_to_fam = [4, 3, 2, 1]

        expected(:, 1, 1) = [70.0_real64, 80.0_real64]
        expected(:, 2, 1) = [-69.0_real64, -78.0_real64]
        expected(:, 1, 2) = [50.0_real64, 60.0_real64]
        expected(:, 2, 2) = [-47.0_real64, -56.0_real64]
        expected(:, 1, 3) = [30.0_real64, 40.0_real64]
        expected(:, 2, 3) = [-25.0_real64, -34.0_real64]
        expected(:, 1, 4) = [10.0_real64, 20.0_real64]
        expected(:, 2, 4) = [-3.0_real64, -12.0_real64]
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_one_family_per_gene")
    end subroutine test_compute_shift_vector_field_one_family_per_gene

    !> Kept from the old suite (`test_zero_distance`), which compared 6 of the 12 elements: a
    !| gene that sits at its family's centroid has the shift exactly 0. Genes [1, 2, 3] and
    !| [4, 5, 6], each in the family with the same centroid.
    subroutine test_compute_shift_vector_field_gene_at_centroid()
        integer(int32), parameter :: n_tissues = 3, n_genes = 2, n_families = 2
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        expression_vectors(:, 2) = [4.0_real64, 5.0_real64, 6.0_real64]
        family_centroids = expression_vectors
        gene_to_fam = [1, 2]

        expected(:, 1, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        expected(:, 2, 1) = 0.0_real64
        expected(:, 1, 2) = [4.0_real64, 5.0_real64, 6.0_real64]
        expected(:, 2, 2) = 0.0_real64
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_gene_at_centroid")
    end subroutine test_compute_shift_vector_field_gene_at_centroid

    !> A gene with the sentinel 0 has no family, so both of its vectors are zero -- not its
    !| expression (a shift from a zero centroid), and nothing read from `family_centroids`.
    !| One family, so a centroid index 0 would read before the array (it did until 307dfd66).
    !| Unassigned genes first, in the middle and last, around two assigned ones:
    !| genes [9, -9], [1, 2], [8, 8], [3, 5], [-7, 6] in [0, 1, 0, 1, 0], family [2, 4]:
    !| - genes 1, 3, 5: centroid [0, 0], shift [0, 0];
    !| - gene 2: centroid [2, 4], shift [-1, -2];
    !| - gene 4: centroid [2, 4], shift [1, 1].
    subroutine test_compute_shift_vector_field_unassigned_genes()
        integer(int32), parameter :: n_tissues = 2, n_genes = 5, n_families = 1
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [9.0_real64, -9.0_real64]
        expression_vectors(:, 2) = [1.0_real64, 2.0_real64]
        expression_vectors(:, 3) = [8.0_real64, 8.0_real64]
        expression_vectors(:, 4) = [3.0_real64, 5.0_real64]
        expression_vectors(:, 5) = [-7.0_real64, 6.0_real64]
        family_centroids(:, 1) = [2.0_real64, 4.0_real64]
        gene_to_fam = [UNASSIGNED, 1_int32, UNASSIGNED, 1_int32, UNASSIGNED]

        expected = 0.0_real64
        expected(:, 1, 2) = [2.0_real64, 4.0_real64]
        expected(:, 2, 2) = [-1.0_real64, -2.0_real64]
        expected(:, 1, 4) = [2.0_real64, 4.0_real64]
        expected(:, 2, 4) = [1.0_real64, 1.0_real64]
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_unassigned_genes")
    end subroutine test_compute_shift_vector_field_unassigned_genes

    !> No gene in any family: every vector is zero, whatever the expression and the centroids.
    !| Genes [1, 2, 3] and [-4, 5, -6], two families.
    subroutine test_compute_shift_vector_field_all_unassigned()
        integer(int32), parameter :: n_tissues = 3, n_genes = 2, n_families = 2
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        expression_vectors(:, 2) = [-4.0_real64, 5.0_real64, -6.0_real64]
        family_centroids(:, 1) = [10.0_real64, 10.0_real64, 10.0_real64]
        family_centroids(:, 2) = [-10.0_real64, -10.0_real64, -10.0_real64]
        gene_to_fam = UNASSIGNED

        expected = 0.0_real64
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_all_unassigned")
    end subroutine test_compute_shift_vector_field_all_unassigned

    !> Four families and two genes, so the centroid matrix is wider than the expression matrix
    !| and two families are left unused. Genes [1, 1] and [2, 2] in families [4, 2]; families
    !| [10, 10], [20, 20], [30, 30], [40, 40]:
    !| - gene 1, family 4: centroid [40, 40], shift [-39, -39];
    !| - gene 2, family 2: centroid [20, 20], shift [-18, -18].
    subroutine test_compute_shift_vector_field_unused_families()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 4
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = 1.0_real64
        expression_vectors(:, 2) = 2.0_real64
        family_centroids(:, 1) = 10.0_real64
        family_centroids(:, 2) = 20.0_real64
        family_centroids(:, 3) = 30.0_real64
        family_centroids(:, 4) = 40.0_real64
        gene_to_fam = [4, 2]

        expected(:, 1, 1) = 40.0_real64
        expected(:, 2, 1) = -39.0_real64
        expected(:, 1, 2) = 20.0_real64
        expected(:, 2, 2) = -18.0_real64
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_unused_families")
    end subroutine test_compute_shift_vector_field_unused_families

    !> 1000 genes over 3 tissues and 7 families, every fifth gene unassigned. The inputs are
    !| built from the result: family f's centroid is 1000*f in every tissue, and gene g in family
    !| f = mod(g, 7) + 1 has the expression 1000*f + [g, 2*g, -g], so its shift is
    !| [g, 2*g, -g] by construction. Every value is an integer below 2**53, so exact.
    subroutine test_compute_shift_vector_field_many_genes()
        integer(int32), parameter :: n_tissues = 3, n_genes = 1000, n_families = 7
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes), i_gene, i_family
        real(real64) :: shift(n_tissues), centroid

        do i_family = 1, n_families
            family_centroids(:, i_family) = 1000.0_real64*i_family
        end do
        do i_gene = 1, n_genes
            shift(1) = real(i_gene, real64)
            shift(2) = 2.0_real64*i_gene
            shift(3) = -real(i_gene, real64)
            if (mod(i_gene, 5) == 0) then
                gene_to_fam(i_gene) = UNASSIGNED
                expression_vectors(:, i_gene) = shift
                expected(:, :, i_gene) = 0.0_real64
            else
                gene_to_fam(i_gene) = mod(i_gene, n_families) + 1
                centroid = 1000.0_real64*gene_to_fam(i_gene)
                expression_vectors(:, i_gene) = centroid + shift
                expected(:, 1, i_gene) = centroid
                expected(:, 2, i_gene) = shift
            end if
        end do
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_many_genes")
    end subroutine test_compute_shift_vector_field_many_genes

    !> The largest finite magnitudes. One gene over four tissues, family 1:
    !| - expression [huge, -huge, 0, huge/2], centroid [huge, 0, -huge, -huge/2];
    !| - shift [huge - huge, -huge - 0, 0 + huge, huge/2 + huge/2] = [0, -huge, huge, huge],
    !|   each exact (huge/2 is exact, and their sum rounds back to huge).
    !| And one difference past `huge`: expression huge, centroid -huge. The true shift 2*huge is
    !| not representable and rounds to +Inf, with `ierr` OK.
    !| QUESTION: pinned as today's behaviour. Finite input gives Inf only where the true value
    !| exceeds `huge`, which the overflow rule of #223 allows; the procedures that read the shift
    !| vectors refuse it with ERR_NAN_INF.
    subroutine test_compute_shift_vector_field_extreme_magnitudes()
        integer(int32), parameter :: n_tissues = 4, n_genes = 1, n_families = 1
        real(real64), parameter :: HUGE_VALUE = huge(1.0_real64)
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: expected(n_tissues, 2, n_genes)
        real(real64) :: overflowing_expression(1, 1), overflowing_centroid(1, 1), overflowing_expected(1, 2, 1)
        integer(int32) :: gene_to_fam(n_genes)

        expression_vectors(:, 1) = [HUGE_VALUE, -HUGE_VALUE, 0.0_real64, HUGE_VALUE/2]
        family_centroids(:, 1) = [HUGE_VALUE, 0.0_real64, -HUGE_VALUE, -HUGE_VALUE/2]
        gene_to_fam = 1
        expected(:, 1, 1) = [HUGE_VALUE, 0.0_real64, -HUGE_VALUE, -HUGE_VALUE/2]
        expected(:, 2, 1) = [0.0_real64, -HUGE_VALUE, HUGE_VALUE, HUGE_VALUE]
        call check_shift_vectors(n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, &
                                 expected, "test_compute_shift_vector_field_extreme_magnitudes: finite")

        overflowing_expression = HUGE_VALUE
        overflowing_centroid = -HUGE_VALUE
        overflowing_expected(1, 1, 1) = -HUGE_VALUE
        overflowing_expected(1, 2, 1) = ieee_value(1.0_real64, ieee_positive_inf)
        call check_shift_vectors(1, 1, 1, overflowing_expression, overflowing_centroid, gene_to_fam, &
                                 overflowing_expected, "test_compute_shift_vector_field_extreme_magnitudes: past huge")
    end subroutine test_compute_shift_vector_field_extreme_magnitudes

    !> The smallest valid sizes, 1 tissue, 1 gene, 1 family: expression [7], centroid [10] gives
    !| centroid [10], shift [-3]; the same gene unassigned gives [0] and [0].
    subroutine test_compute_shift_vector_field_smallest_sizes()
        real(real64) :: expression_vectors(1, 1), family_centroids(1, 1), expected(1, 2, 1)
        integer(int32) :: gene_to_fam(1)

        expression_vectors = 7.0_real64
        family_centroids = 10.0_real64

        gene_to_fam = 1
        expected(1, :, 1) = [10.0_real64, -3.0_real64]
        call check_shift_vectors(1, 1, 1, expression_vectors, family_centroids, gene_to_fam, expected, &
                                 "test_compute_shift_vector_field_smallest_sizes: in the family")

        gene_to_fam = UNASSIGNED
        expected = 0.0_real64
        call check_shift_vectors(1, 1, 1, expression_vectors, family_centroids, gene_to_fam, expected, &
                                 "test_compute_shift_vector_field_smallest_sizes: unassigned")
    end subroutine test_compute_shift_vector_field_smallest_sizes

    !> `gene_to_fam` (argument 6) holds ids in [1, n_families] or the sentinel 0; 1, n_families
    !| and 0 are accepted in the cases above. Kept from the old suite
    !| (`test_invalid_family_mapping`): n_families + 1 is ERR_INVALID_INPUT; so are -1 and
    !| huge. Each sits in the last gene, so the check must cover every gene.
    subroutine test_compute_shift_vector_field_family_id_range()
        integer(int32), parameter :: n_tissues = 2, n_genes = 3, n_families = 3
        integer(int32), parameter :: N_OUT_OF_RANGE = 3
        integer(int32) :: out_of_range(N_OUT_OF_RANGE)
        character(len=14) :: out_of_range_names(N_OUT_OF_RANGE)
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: shift_vectors(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_id

        out_of_range(1) = n_families + 1_int32
        out_of_range(2) = -1_int32
        out_of_range(3) = huge(1_int32)
        out_of_range_names(1) = "n_families + 1"
        out_of_range_names(2) = "-1"
        out_of_range_names(3) = "huge"
        expression_vectors = 1.0_real64
        family_centroids = 2.0_real64
        do i_id = 1, N_OUT_OF_RANGE
            gene_to_fam(1) = 1_int32
            gene_to_fam(2) = n_families
            gene_to_fam(3) = out_of_range(i_id)
            call compute_shift_vector_field(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                            gene_to_fam, shift_vectors, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_shift_vector_field_family_id_range: " &
                            //trim(out_of_range_names(i_id)), 6)
        end do
    end subroutine test_compute_shift_vector_field_family_id_range

    !> Each extent is ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1, naming it:
    !| `n_tissues` (argument 1), `n_genes` (argument 2; 0 kept from the old suite's
    !| `test_dimension_edge_cases`) and `n_families` (argument 3). 1 is in
    !| `test_compute_shift_vector_field_smallest_sizes`.
    subroutine test_compute_shift_vector_field_dimensions()
        real(real64) :: expression_vectors(1, 1), family_centroids(1, 1), shift_vectors(1, 2, 1)
        integer(int32) :: gene_to_fam(1), ierr

        expression_vectors = 1.0_real64
        family_centroids = 1.0_real64
        gene_to_fam = 1

        call compute_shift_vector_field(0, 1, 1, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_shift_vector_field_dimensions: n_tissues 0", 1)
        call compute_shift_vector_field(-1, 1, 1, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_shift_vector_field_dimensions: n_tissues -1", 1)

        call compute_shift_vector_field(1, 0, 1, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_shift_vector_field_dimensions: n_genes 0", 2)
        call compute_shift_vector_field(1, -1, 1, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_shift_vector_field_dimensions: n_genes -1", 2)

        ! the gene's id 1 is past n_families as well; the extent is checked first and named
        call compute_shift_vector_field(1, 1, 0, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_shift_vector_field_dimensions: n_families 0", 3)
        call compute_shift_vector_field(1, 1, -1, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_shift_vector_field_dimensions: n_families -1", 3)
    end subroutine test_compute_shift_vector_field_dimensions

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `expression_vectors` (argument 4) and in
    !| `family_centroids` (argument 5). Each sits in the last element, so the check must cover
    !| every element; in `expression_vectors` it belongs to an unassigned gene, whose expression
    !| the computation never reads.
    subroutine test_compute_shift_vector_field_rejects_nan_and_inf()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 2
        real(real64) :: expression_vectors(n_tissues, n_genes), family_centroids(n_tissues, n_families)
        real(real64) :: shift_vectors(n_tissues, 2, n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_value

        gene_to_fam = [2_int32, UNASSIGNED]
        do i_value = 1, N_NON_FINITE
            expression_vectors = 1.0_real64
            family_centroids = 2.0_real64
            expression_vectors(n_tissues, n_genes) = non_finite_value(i_value)
            call compute_shift_vector_field(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                            gene_to_fam, shift_vectors, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_shift_vector_field_rejects_nan_and_inf: expression_vectors " &
                            //trim(NON_FINITE_NAMES(i_value)), 4)

            expression_vectors = 1.0_real64
            family_centroids(n_tissues, n_families) = non_finite_value(i_value)
            call compute_shift_vector_field(n_tissues, n_genes, n_families, expression_vectors, family_centroids, &
                                            gene_to_fam, shift_vectors, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_shift_vector_field_rejects_nan_and_inf: family_centroids " &
                            //trim(NON_FINITE_NAMES(i_value)), 5)
        end do
    end subroutine test_compute_shift_vector_field_rejects_nan_and_inf

end module mod_test_compute_shift_vector_field
