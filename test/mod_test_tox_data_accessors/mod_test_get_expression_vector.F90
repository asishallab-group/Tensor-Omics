!> The `get_expression_vector` cases: gene j's vector is column j of `expression_vectors`; a gene
!| index outside 1..n_genes is refused with ERR_INVALID_INPUT, an `out_vector` of another length
!| than the samples with ERR_SIZE_MISMATCH, in that order; a refusal leaves `out_vector` alone.
module mod_test_get_expression_vector
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use tox_data_accessors, only: get_expression_vector
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: N_SAMPLES, N_GENES, UNWRITTEN, expression_vectors
    implicit none
    private
    public :: get_all_tests_get_expression_vector

contains

    !> Every case of `get_expression_vector`.
    function get_all_tests_get_expression_vector() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(5))
        all_tests(1) = test_case("test_get_expression_vector_values", test_get_expression_vector_values)
        all_tests(2) = test_case("test_get_expression_vector_gene_bound", test_get_expression_vector_gene_bound)
        all_tests(3) = test_case("test_get_expression_vector_length", test_get_expression_vector_length)
        all_tests(4) = test_case("test_get_expression_vector_no_sample", test_get_expression_vector_no_sample)
        all_tests(5) = test_case("test_get_expression_vector_copies_non_finite", &
                                 test_get_expression_vector_copies_non_finite)
    end function get_all_tests_get_expression_vector

    !> Each gene's column, [10*j + 1, 10*j + 2]: [11, 12], [21, 22], [31, 32].
    subroutine test_get_expression_vector_values()
        real(real64) :: expression(N_SAMPLES, N_GENES), vector(N_SAMPLES), expected(N_SAMPLES, N_GENES)
        integer(int32) :: i_gene, ierr

        expression = expression_vectors()
        expected = reshape([11.0_real64, 12.0_real64, 21.0_real64, 22.0_real64, 31.0_real64, 32.0_real64], [N_SAMPLES, N_GENES])
        do i_gene = 1, N_GENES
            vector = UNWRITTEN
            call get_expression_vector(i_gene, expression, vector, ierr)
            call assert_err(ierr, ERR_OK, "test_get_expression_vector_values: ierr of gene "//i_gene)
            call assert_equal_array_real(vector, expected(:, i_gene), N_SAMPLES, 0.0_real64, &
                                         "test_get_expression_vector_values: gene "//i_gene)
        end do
    end subroutine test_get_expression_vector_values

    !> 0 and n_genes + 1 are refused with ERR_INVALID_INPUT and leave `out_vector` alone; so
    !| is the gene index when the vector's length is wrong too, which is checked second.
    subroutine test_get_expression_vector_gene_bound()
        real(real64) :: expression(N_SAMPLES, N_GENES), vector(N_SAMPLES), expected(N_SAMPLES), too_long(N_SAMPLES + 1)
        integer(int32) :: ierr

        expression = expression_vectors()
        expected = UNWRITTEN
        vector = UNWRITTEN
        call get_expression_vector(0, expression, vector, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_expression_vector_gene_bound: 0")
        call assert_equal_array_real(vector, expected, N_SAMPLES, 0.0_real64, "test_get_expression_vector_gene_bound: 0 writes")
        call get_expression_vector(N_GENES + 1, expression, vector, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_expression_vector_gene_bound: n_genes + 1")
        call assert_equal_array_real(vector, expected, N_SAMPLES, 0.0_real64, &
                                     "test_get_expression_vector_gene_bound: n_genes + 1 writes")
        call get_expression_vector(0, expression, too_long, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_expression_vector_gene_bound: 0, and a vector too long")
    end subroutine test_get_expression_vector_gene_bound

    !> An `out_vector` one longer or one shorter than the samples is refused with
    !| ERR_SIZE_MISMATCH.
    subroutine test_get_expression_vector_length()
        real(real64) :: expression(N_SAMPLES, N_GENES), too_long(N_SAMPLES + 1), too_short(N_SAMPLES - 1)
        integer(int32) :: ierr

        expression = expression_vectors()
        call get_expression_vector(1, expression, too_long, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_expression_vector_length: one too long")
        call get_expression_vector(1, expression, too_short, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_expression_vector_length: one too short")
    end subroutine test_get_expression_vector_length

    !> Genes with no sample have a vector of no value: an empty `out_vector` is the right length.
    subroutine test_get_expression_vector_no_sample()
        real(real64) :: expression(0, N_GENES), vector(0)
        integer(int32) :: ierr

        call get_expression_vector(1, expression, vector, ierr)
        call assert_err(ierr, ERR_OK, "test_get_expression_vector_no_sample")
    end subroutine test_get_expression_vector_no_sample

    !> The accessors copy what is stored as it is, NaN and Inf included, with ERR_OK, as the
    !| module documents: they check indices and shapes, not values.
    subroutine test_get_expression_vector_copies_non_finite()
        real(real64) :: expression(N_SAMPLES, N_GENES), vector(N_SAMPLES)
        integer(int32) :: ierr

        expression = expression_vectors()
        expression(1, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
        expression(2, 2) = ieee_value(1.0_real64, ieee_positive_inf)
        call get_expression_vector(2, expression, vector, ierr)
        call assert_err(ierr, ERR_OK, "test_get_expression_vector_copies_non_finite: ierr")
        call assert_true(ieee_is_nan(vector(1)), "test_get_expression_vector_copies_non_finite: the NaN")
        call assert_true(vector(2) == ieee_value(1.0_real64, ieee_positive_inf), &
                         "test_get_expression_vector_copies_non_finite: the +Inf")
    end subroutine test_get_expression_vector_copies_non_finite
end module mod_test_get_expression_vector
