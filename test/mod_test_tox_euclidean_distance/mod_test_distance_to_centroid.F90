!> The `distance_to_centroid` cases: for every gene, the Euclidean distance to its own family's
!| centroid; genes in no order of their families, a gene at its centroid, unassigned genes (the
!| sentinel 0, which has no centroid and gets the distance sentinel -1), more families than genes,
!| many tissues, many genes, the largest finite magnitudes, the smallest sizes, the family-id
!| range, the dimension checks and the NaN/Inf checks.
!|
!| Every expected value is written out by hand. Each gene is its centroid plus an offset that is
!| a Pythagorean triple (or the zero vector), so every distance is an integer, exact, and every
!| comparison has tolerance 0. The offsets differ from tissue to tissue, so that reading the
!| matrices transposed, or the wrong family's centroid, gives another number.
module mod_test_distance_to_centroid
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use tox_euclidean_distance, only: distance_to_centroid
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_euclidean_distance_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_distance_to_centroid

    !> `gene_to_fam` of a gene in no family (M_GENE_TO_FAM_SENTINEL in src/macros.h).
    integer(int32), parameter :: UNASSIGNED = 0_int32
    !> The distance such a gene gets (CM_DISTANCE_SENTINEL in tox_euclidean_distance_impl).
    real(real64), parameter :: NO_DISTANCE = -1.0_real64

contains

    !> Every case of `distance_to_centroid`.
    function get_all_tests_distance_to_centroid() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(13))
        all_tests(1) = test_case("test_distance_to_centroid_values", test_distance_to_centroid_values)
        all_tests(2) = test_case("test_distance_to_centroid_smallest_sizes", test_distance_to_centroid_smallest_sizes)
        all_tests(3) = test_case("test_distance_to_centroid_unassigned_genes", test_distance_to_centroid_unassigned_genes)
        all_tests(4) = test_case("test_distance_to_centroid_all_unassigned", test_distance_to_centroid_all_unassigned)
        all_tests(5) = test_case("test_distance_to_centroid_unused_families", test_distance_to_centroid_unused_families)
        all_tests(6) = test_case("test_distance_to_centroid_many_tissues", test_distance_to_centroid_many_tissues)
        all_tests(7) = test_case("test_distance_to_centroid_many_genes", test_distance_to_centroid_many_genes)
        all_tests(8) = test_case("test_distance_to_centroid_largest_magnitudes", test_distance_to_centroid_largest_magnitudes)
        all_tests(9) = test_case("test_distance_to_centroid_family_id_range", test_distance_to_centroid_family_id_range)
        all_tests(10) = test_case("test_distance_to_centroid_dimensions", test_distance_to_centroid_dimensions)
        all_tests(11) = test_case("test_distance_to_centroid_rejects_nan_and_inf", &
                                  test_distance_to_centroid_rejects_nan_and_inf)
        all_tests(12) = test_case("test_distance_to_centroid_squares_past_huge", test_distance_to_centroid_squares_past_huge)
        all_tests(13) = test_case("test_distance_to_centroid_squares_below_tiny", test_distance_to_centroid_squares_below_tiny)
    end function get_all_tests_distance_to_centroid

    !> Calls the procedure and checks `ierr` OK and every element of `distances` exactly.
    !| `distances` starts at UNWRITTEN and has one gene more than the call is told, so an element
    !| the procedure skips and a write past the last gene both show.
    subroutine check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, msg)
        integer(int32), intent(in) :: n_tissues, n_genes, n_families
        real(real64), intent(in) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
        real(real64), intent(in) :: expected(n_genes)
        character(*), intent(in) :: msg
        real(real64) :: distances(n_genes + 1)
        integer(int32) :: ierr

        distances = UNWRITTEN
        call distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, n_tissues, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(distances, expected, n_genes, 0.0_real64, msg//": distances")
        call assert_equal_real(distances(n_genes + 1), UNWRITTEN, 0.0_real64, msg//": no write past the last gene")
    end subroutine check_distances

    !> Kept from the old suite (`test_distance_to_centroid_basic` and
    !| `test_column_major_optimization`, whose expected values sqrt(0.5), sqrt(4.5) and
    !| sqrt(252525.25) are not exact): 3 tissues, 5 genes, 3 families, the genes in no order of
    !| their families, every family used, the first and the last id. Families 1..3 are
    !| [1, 2, 3], [-4, 0, 4], [10, 20, 30]; gene_to_fam = [2, 3, 1, 3, 1]:
    !| - gene 1, family 2, offset [1, 2, 2]: gene [-3, 2, 6], distance 3;
    !| - gene 2, family 3, offset [2, 3, 6]: gene [12, 23, 36], distance 7;
    !| - gene 3, family 1, offset [-1, -4, 8]: gene [0, -2, 11], distance 9;
    !| - gene 4, family 3, at its centroid: gene [10, 20, 30], distance 0;
    !| - gene 5, family 1, offset [2, -10, 11]: gene [3, -8, 14], distance 15.
    subroutine test_distance_to_centroid_values()
        integer(int32), parameter :: n_tissues = 3, n_genes = 5, n_families = 3
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        centroids(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        centroids(:, 2) = [-4.0_real64, 0.0_real64, 4.0_real64]
        centroids(:, 3) = [10.0_real64, 20.0_real64, 30.0_real64]
        genes(:, 1) = [-3.0_real64, 2.0_real64, 6.0_real64]
        genes(:, 2) = [12.0_real64, 23.0_real64, 36.0_real64]
        genes(:, 3) = [0.0_real64, -2.0_real64, 11.0_real64]
        genes(:, 4) = [10.0_real64, 20.0_real64, 30.0_real64]
        genes(:, 5) = [3.0_real64, -8.0_real64, 14.0_real64]
        gene_to_fam = [2, 3, 1, 3, 1]

        expected = [3.0_real64, 7.0_real64, 9.0_real64, 0.0_real64, 15.0_real64]
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_values")
    end subroutine test_distance_to_centroid_values

    !> Kept from the old suite (`test_distance_to_centroid_single`, 2 tissues, a gene at its
    !| centroid, tolerance 1e-15): one gene in one family. [1, 2] at the centroid [1, 2] is at
    !| distance 0; the smallest valid sizes, 1 tissue, 1 gene, 1 family, [7] to the centroid [10]
    !| is 3.
    subroutine test_distance_to_centroid_smallest_sizes()
        real(real64) :: genes(2, 1), centroids(2, 1), single_gene(1, 1), single_centroid(1, 1), expected(1)
        integer(int32) :: gene_to_fam(1)

        gene_to_fam = 1
        genes(:, 1) = [1.0_real64, 2.0_real64]
        centroids = genes
        expected = 0.0_real64
        call check_distances(2, 1, 1, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_smallest_sizes: one gene at its centroid")

        single_gene = 7.0_real64
        single_centroid = 10.0_real64
        expected = 3.0_real64
        call check_distances(1, 1, 1, single_gene, single_centroid, gene_to_fam, expected, &
                             "test_distance_to_centroid_smallest_sizes: one tissue")
    end subroutine test_distance_to_centroid_smallest_sizes

    !> Kept from the old suite (`test_distance_to_centroid_invalid_families`, second call): a gene
    !| with the sentinel 0 has no centroid and gets the distance sentinel -1, which no distance can
    !| be. One family, so a centroid index 0 would read before the array. Unassigned genes first,
    !| in the middle and last, around two assigned ones: genes [9, -9], [4, 7], [8, 8], [-11, 8],
    !| [-7, 6] in [0, 1, 0, 1, 0], family [1, 3]:
    !| - genes 1, 3, 5: -1;
    !| - gene 2, offset [3, 4]: distance 5;
    !| - gene 4, offset [-12, 5]: distance 13.
    subroutine test_distance_to_centroid_unassigned_genes()
        integer(int32), parameter :: n_tissues = 2, n_genes = 5, n_families = 1
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        centroids(:, 1) = [1.0_real64, 3.0_real64]
        genes(:, 1) = [9.0_real64, -9.0_real64]
        genes(:, 2) = [4.0_real64, 7.0_real64]
        genes(:, 3) = [8.0_real64, 8.0_real64]
        genes(:, 4) = [-11.0_real64, 8.0_real64]
        genes(:, 5) = [-7.0_real64, 6.0_real64]
        gene_to_fam = [UNASSIGNED, 1_int32, UNASSIGNED, 1_int32, UNASSIGNED]

        expected = [NO_DISTANCE, 5.0_real64, NO_DISTANCE, 13.0_real64, NO_DISTANCE]
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_unassigned_genes")
    end subroutine test_distance_to_centroid_unassigned_genes

    !> No gene in any family: every distance is the sentinel -1, whatever the expression and the
    !| centroids. Genes [1, 2, 3] and [-4, 5, -6], two families.
    subroutine test_distance_to_centroid_all_unassigned()
        integer(int32), parameter :: n_tissues = 3, n_genes = 2, n_families = 2
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        genes(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        genes(:, 2) = [-4.0_real64, 5.0_real64, -6.0_real64]
        centroids(:, 1) = 10.0_real64
        centroids(:, 2) = -10.0_real64
        gene_to_fam = UNASSIGNED

        expected = NO_DISTANCE
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_all_unassigned")
    end subroutine test_distance_to_centroid_all_unassigned

    !> Four families and two genes, so the centroid matrix is wider than the expression matrix and
    !| two families are left unused. Families [10, 10], [20, 20], [30, 30], [40, 40]; genes in
    !| families [4, 2]:
    !| - gene 1, family 4, offset [6, 8]: gene [46, 48], distance 10;
    !| - gene 2, family 2, offset [-5, 12]: gene [15, 32], distance 13.
    subroutine test_distance_to_centroid_unused_families()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 4
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        centroids(:, 1) = 10.0_real64
        centroids(:, 2) = 20.0_real64
        centroids(:, 3) = 30.0_real64
        centroids(:, 4) = 40.0_real64
        genes(:, 1) = [46.0_real64, 48.0_real64]
        genes(:, 2) = [15.0_real64, 32.0_real64]
        gene_to_fam = [4, 2]

        expected = [10.0_real64, 13.0_real64]
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_unused_families")
    end subroutine test_distance_to_centroid_unused_families

    !> Kept from the old suite (`test_distance_to_centroid_high_dim`): 100 tissues, 3 genes, one
    !| family. The centroid is [2, 3, ..., 101]; genes 1, 2, 3 are the centroid minus 1, the
    !| centroid itself and the centroid plus 1 in every tissue, so their sums of squares are 100,
    !| 0 and 100, and the distances exactly 10, 0 and 10.
    subroutine test_distance_to_centroid_many_tissues()
        integer(int32), parameter :: n_tissues = 100, n_genes = 3, n_families = 1
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes), i_tissue

        do i_tissue = 1, n_tissues
            genes(i_tissue, 1) = real(i_tissue, real64)
            genes(i_tissue, 2) = real(i_tissue + 1, real64)
            genes(i_tissue, 3) = real(i_tissue + 2, real64)
            centroids(i_tissue, 1) = real(i_tissue + 1, real64)
        end do
        gene_to_fam = 1

        expected = [10.0_real64, 0.0_real64, 10.0_real64]
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_many_tissues")
    end subroutine test_distance_to_centroid_many_tissues

    !> In place of the old `test_performance_large_scale` (1000 genes of trigonometric values,
    !| which asserted only that some distance is positive): 1000 genes over 3 tissues and 7
    !| families, every fifth gene unassigned. The inputs are built from the result: family f's
    !| centroid is [100*f, -100*f, 1000*f], and gene g in family f = mod(g, 7) + 1 is that centroid
    !| plus [g, 2*g, 2*g], a triple times g, so its distance is 3*g; an unassigned gene gets -1.
    !| Every value is an integer below 2**53, so exact.
    subroutine test_distance_to_centroid_many_genes()
        integer(int32), parameter :: n_tissues = 3, n_genes = 1000, n_families = 7
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes), i_gene, i_family
        real(real64) :: offset(n_tissues)

        do i_family = 1, n_families
            centroids(:, i_family) = [100.0_real64, -100.0_real64, 1000.0_real64]*i_family
        end do
        do i_gene = 1, n_genes
            offset = [1.0_real64, 2.0_real64, 2.0_real64]*i_gene
            if (mod(i_gene, 5) == 0) then
                gene_to_fam(i_gene) = UNASSIGNED
                genes(:, i_gene) = offset
                expected(i_gene) = NO_DISTANCE
            else
                gene_to_fam(i_gene) = mod(i_gene, n_families) + 1
                genes(:, i_gene) = centroids(:, gene_to_fam(i_gene)) + offset
                expected(i_gene) = 3.0_real64*i_gene
            end if
        end do
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_many_genes")
    end subroutine test_distance_to_centroid_many_genes

    !> The largest finite magnitudes whose differences need no square past `huge`, two tissues:
    !| - gene 1 [huge, -huge] at its centroid [huge, -huge]: distance 0;
    !| - gene 2 [huge, 0] in a family at [-huge, 0]: the true distance 2*huge is past `huge`, and
    !|   the result is +Inf with `ierr` OK, as for `euclidean_distance` (see
    !|   `test_euclidean_distance_largest_magnitudes`).
    subroutine test_distance_to_centroid_largest_magnitudes()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 2
        real(real64), parameter :: HUGE_VALUE = huge(1.0_real64)
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        genes(:, 1) = [HUGE_VALUE, -HUGE_VALUE]
        genes(:, 2) = [HUGE_VALUE, 0.0_real64]
        centroids(:, 1) = [HUGE_VALUE, -HUGE_VALUE]
        centroids(:, 2) = [-HUGE_VALUE, 0.0_real64]
        gene_to_fam = [1, 2]

        expected(1) = 0.0_real64
        expected(2) = ieee_value(1.0_real64, ieee_positive_inf)
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_largest_magnitudes")
    end subroutine test_distance_to_centroid_largest_magnitudes

    !> `gene_to_fam` (argument 5) holds ids in [1, n_families] or the sentinel 0; 1, n_families
    !| and 0 are accepted in the cases above. Kept from the old suite
    !| (`test_distance_to_centroid_invalid_families`, first call, and the Python and R suites'
    !| negative ids): n_families + 1 is ERR_INVALID_INPUT; so are -1, huge and -huge. Each sits in
    !| the last gene, so the check must cover every gene.
    subroutine test_distance_to_centroid_family_id_range()
        integer(int32), parameter :: n_tissues = 2, n_genes = 3, n_families = 2
        integer(int32), parameter :: N_OUT_OF_RANGE = 4
        integer(int32) :: out_of_range(N_OUT_OF_RANGE)
        character(len=14) :: out_of_range_names(N_OUT_OF_RANGE)
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), distances(n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_id

        out_of_range = [n_families + 1_int32, -1_int32, huge(1_int32), -huge(1_int32)]
        out_of_range_names(1) = "n_families + 1"
        out_of_range_names(2) = "-1"
        out_of_range_names(3) = "huge"
        out_of_range_names(4) = "-huge"
        genes = 1.0_real64
        centroids = 2.0_real64
        do i_id = 1, N_OUT_OF_RANGE
            gene_to_fam(1) = 1_int32
            gene_to_fam(2) = n_families
            gene_to_fam(3) = out_of_range(i_id)
            call distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, n_tissues, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_distance_to_centroid_family_id_range: " &
                            //trim(out_of_range_names(i_id)), 5)
        end do
    end subroutine test_distance_to_centroid_family_id_range

    !> Each extent is ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1, naming it: `n_genes`
    !| (argument 1; 0 kept from the old Python and R suites' empty input), `n_families`
    !| (argument 2) and `n_tissues` (argument 7). 1 is in
    !| `test_distance_to_centroid_smallest_sizes`.
    subroutine test_distance_to_centroid_dimensions()
        real(real64) :: genes(1, 1), centroids(1, 1), distances(1)
        integer(int32) :: gene_to_fam(1), ierr

        genes = 1.0_real64
        centroids = 1.0_real64
        gene_to_fam = 1

        call distance_to_centroid(0, 1, genes, centroids, gene_to_fam, distances, 1, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_distance_to_centroid_dimensions: n_genes 0", 1)
        call distance_to_centroid(-1, 1, genes, centroids, gene_to_fam, distances, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_distance_to_centroid_dimensions: n_genes -1", 1)

        ! the gene's id 1 is past n_families as well; the extent is checked first and named
        call distance_to_centroid(1, 0, genes, centroids, gene_to_fam, distances, 1, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_distance_to_centroid_dimensions: n_families 0", 2)
        call distance_to_centroid(1, -1, genes, centroids, gene_to_fam, distances, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_distance_to_centroid_dimensions: n_families -1", 2)

        call distance_to_centroid(1, 1, genes, centroids, gene_to_fam, distances, 0, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_distance_to_centroid_dimensions: n_tissues 0", 7)
        call distance_to_centroid(1, 1, genes, centroids, gene_to_fam, distances, -1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_distance_to_centroid_dimensions: n_tissues -1", 7)
    end subroutine test_distance_to_centroid_dimensions

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `genes` (argument 3) and in `centroids`
    !| (argument 4). Each sits in the last element, so the check must cover every element, and
    !| belongs to what the computation never reads: an unassigned gene, a family no gene is in.
    subroutine test_distance_to_centroid_rejects_nan_and_inf()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 2
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), distances(n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_value

        gene_to_fam = [1_int32, UNASSIGNED]
        do i_value = 1, N_NON_FINITE
            genes = 1.0_real64
            centroids = 2.0_real64
            genes(n_tissues, n_genes) = non_finite_value(i_value)
            call distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, n_tissues, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_distance_to_centroid_rejects_nan_and_inf: genes " &
                            //trim(NON_FINITE_NAMES(i_value)), 3)

            genes = 1.0_real64
            centroids(n_tissues, n_families) = non_finite_value(i_value)
            call distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, n_tissues, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_distance_to_centroid_rejects_nan_and_inf: centroids " &
                            //trim(NON_FINITE_NAMES(i_value)), 4)
        end do
    end subroutine test_distance_to_centroid_rejects_nan_and_inf

    !> BUG: distances whose squares overflow are +Inf, with `ierr` OK, although they are finite
    !| (see `test_euclidean_distance_squares_past_huge`: the procedure takes the same unscaled
    !| sum of squares). Two tissues, two genes, two families:
    !| - gene 1 [4*2**600, 3*2**600] in family [2**600, -2**600], offset [3, 4]*2**600: distance
    !|   5*2**600; today +Inf;
    !| - gene 2 [huge, 0] in family [0, 0]: distance huge; today +Inf.
    subroutine test_distance_to_centroid_squares_past_huge()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 2
        real(real64), parameter :: LARGE_SCALE = 2.0_real64**600
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        centroids(:, 1) = [1.0_real64, -1.0_real64]*LARGE_SCALE
        centroids(:, 2) = 0.0_real64
        genes(:, 1) = [4.0_real64, 3.0_real64]*LARGE_SCALE
        genes(:, 2) = [huge(1.0_real64), 0.0_real64]
        gene_to_fam = [1, 2]

        expected(1) = 5.0_real64*LARGE_SCALE
        expected(2) = huge(1.0_real64)
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_squares_past_huge")
    end subroutine test_distance_to_centroid_squares_past_huge

    !> BUG: distances whose squares underflow are 0, or wrong, with `ierr` OK (see
    !| `test_euclidean_distance_squares_below_tiny`). Two tissues, two genes, one family
    !| [2**-540, 0]:
    !| - gene 1 [4*2**-540, 4*2**-540], offset [3, 4]*2**-540: distance 5*2**-540; today 0;
    !| - gene 2 [2**-540, 0], at its centroid: distance 0, today too.
    subroutine test_distance_to_centroid_squares_below_tiny()
        integer(int32), parameter :: n_tissues = 2, n_genes = 2, n_families = 1
        real(real64), parameter :: SMALL_SCALE = 2.0_real64**(-540)
        real(real64) :: genes(n_tissues, n_genes), centroids(n_tissues, n_families), expected(n_genes)
        integer(int32) :: gene_to_fam(n_genes)

        centroids(:, 1) = [1.0_real64, 0.0_real64]*SMALL_SCALE
        genes(:, 1) = [4.0_real64, 4.0_real64]*SMALL_SCALE
        genes(:, 2) = centroids(:, 1)
        gene_to_fam = 1

        expected(1) = 5.0_real64*SMALL_SCALE
        expected(2) = 0.0_real64
        call check_distances(n_tissues, n_genes, n_families, genes, centroids, gene_to_fam, expected, &
                             "test_distance_to_centroid_squares_below_tiny")
    end subroutine test_distance_to_centroid_squares_below_tiny

end module mod_test_distance_to_centroid
