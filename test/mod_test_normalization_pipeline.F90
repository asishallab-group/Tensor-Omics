!> Unit test suite for normalization_pipeline routine.
module mod_test_normalization_pipeline
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public


contains

  !> Get array of all available tests.
  function get_all_tests_normalization_pipeline() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(4))
    all_tests(1) = test_case("test_pipeline_basic", test_pipeline_basic)
    all_tests(2) = test_case("test_pipeline_edge_cases", test_pipeline_edge_cases)
    all_tests(3) = test_case("test_pipeline_vs_manual", test_pipeline_vs_manual)
    all_tests(4) = test_case("test_pipeline_empty_matrix", test_pipeline_empty_matrix)
  end function get_all_tests_normalization_pipeline
 
  

  !> Basic test: pipeline with small matrix, no fold change
  subroutine test_pipeline_basic()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_groups = 2
    real(real64), dimension(n_tissues, n_genes) :: expr
    real(real64), dimension(n_groups, n_genes) :: log_transformed_expr
    integer(int32), dimension(n_groups) :: group_sizes
    real(real64) :: span
    integer(int32) :: degree, ierr, i, j

    do i = 1, n_genes
      do j = 1, n_tissues
        expr(j, i) = dble(i) + dble(j) * 0.25d0
      end do
    end do
    group_sizes = [3,3]
    span = 0.75d0
    degree = 2

    call normalization_pipeline_alloc(n_genes, n_tissues, expr, log_transformed_expr, group_sizes, n_groups, span, degree, use_quantile=.true., ierr=ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_basic: normalization_pipeline returned error")
    call assert_no_nan_real(log_transformed_expr, n_genes * n_groups, "test_pipeline_basic: NaN in output")
    call assert_equal_int(size(log_transformed_expr), n_genes * n_groups, "test_pipeline_basic: output size incorrect")
    call assert_true(all(log_transformed_expr >= 0.0d0), "test_pipeline_basic: output should be non-negative")
  end subroutine test_pipeline_basic

  !> Edge case test: zero input, check output shape and values
  subroutine test_pipeline_edge_cases()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_groups = 2
    real(real64), dimension(n_tissues, n_genes) :: expr
    real(real64), dimension(n_groups, n_genes) :: log_transformed_expr
    integer(int32), dimension(n_groups) :: group_sizes
    real(real64) :: span
    integer(int32) :: degree, ierr

    expr = 0.0d0
    group_sizes = [3,3]
    span = 0.75d0
    degree = 2

    call normalization_pipeline_alloc(n_genes, n_tissues, expr, log_transformed_expr, group_sizes, n_groups, span, degree, use_quantile=.true., ierr=ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_pipeline_edge_cases: expected error for zero-variance input")
  end subroutine test_pipeline_edge_cases

  !> Test: pipeline output matches manual stepwise normalization (no fold change)
  subroutine test_pipeline_vs_manual()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_groups = 2
    real(real64), dimension(n_tissues, n_genes) :: expr, stddev, quantile
    real(real64), dimension(n_groups, n_genes) :: log_transformed_expr
    real(real64), dimension(n_groups, n_genes) :: tiss_avg, log2trans, buf_avg_no_quant, log_transformed_expr_no_quant, log2trans_no_quant
    integer(int32), dimension(n_groups) :: group_sizes
    real(real64), dimension(n_genes) :: rank_means, tmp_col
    integer(int32), dimension(n_genes) :: tmp_perm
    real(real64) :: span
    integer(int32) :: degree, ierr, i, j

    ! Build varied data with non-zero variance per gene
    do i = 1, n_genes
      do j = 1, n_tissues
        expr(j, i) = dble(i) * 2.0d0 + dble(j) * 0.5d0
      end do
    end do
    group_sizes = [3,3]
    span = 0.75d0
    degree = 2

    call normalize_by_std_dev_alloc(n_genes, n_tissues, expr, stddev, span, degree, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: normalize_by_std_dev returned error")

    call quantile_normalization(n_genes, n_tissues, stddev, quantile, rank_means, tmp_col, tmp_perm, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: quantile_normalization returned error")
    call calc_tiss_avg(n_genes, n_groups, group_sizes, quantile, tiss_avg, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: calc_tiss_avg returned error")
    call log2_transformation(n_genes, n_groups, tiss_avg, log2trans, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: log2_transformation returned error")

    ! Manual stepwise normalization without quantile
    call calc_tiss_avg(n_genes, n_groups, group_sizes, stddev, buf_avg_no_quant, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: calc_tiss_avg (no quantile) returned error")
    call log2_transformation(n_genes, n_groups, buf_avg_no_quant, log2trans_no_quant, ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: log2_transformation (no quantile) returned error")

    ! Pipeline normalization
    call normalization_pipeline_alloc(n_genes, n_tissues, expr, log_transformed_expr, group_sizes, n_groups, span, degree, use_quantile=.true., ierr=ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: normalization_pipeline returned error")
    call assert_equal_array_real(log_transformed_expr, log2trans, size(log_transformed_expr, kind=int32), 1d-12, "test_pipeline_vs_manual: test_pipeline_vs_manual: pipeline and manual outputs differ")

    call normalization_pipeline_alloc(n_genes, n_tissues, expr, log_transformed_expr_no_quant, group_sizes, n_groups, span, degree, ierr=ierr)
    call assert_equal_int(ierr, ERR_OK, "test_pipeline_vs_manual: normalization_pipeline (no quantile) returned error")
    call assert_equal_array_real(log_transformed_expr_no_quant, log2trans_no_quant, size(log_transformed_expr_no_quant, kind=int32), 1d-12, "test_pipeline_vs_manual: pipeline and manual outputs differ (no quantile)")
  end subroutine test_pipeline_vs_manual

  !> Test with empty input matrix.
  subroutine test_pipeline_empty_matrix()
    integer(int32) :: n_genes, n_tissues, n_groups, degree, ierr
    real(real64) :: span
    real(real64), dimension(1) :: expr, log_transformed_expr
    integer(int32), dimension(1) :: group_sizes
    n_genes = 0; n_tissues = 0; n_groups = 0
    span = 0.75d0
    degree = 2
    call normalization_pipeline_alloc(n_genes, n_tissues, expr, log_transformed_expr, group_sizes, n_groups, span, degree, ierr=ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "normalization_pipeline should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_pipeline_empty_matrix

end module mod_test_normalization_pipeline
