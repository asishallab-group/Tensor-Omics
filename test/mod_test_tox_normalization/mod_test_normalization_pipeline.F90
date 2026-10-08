!> The `normalization_pipeline` cases: the pipeline must equal its steps run one by one, for one,
!| two and three tissues (the pipeline reuses different columns of its output as scratch for each),
!| with and without quantile normalization; then each early exit, and the input checks.
module mod_test_normalization_pipeline
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use tox_normalization
  use mod_test_tox_normalization_fixtures, only: fill_linear_trend
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_normalization_pipeline() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(7))
    all_tests(1) = test_case("test_pipeline_matches_its_steps", test_pipeline_matches_its_steps)
    all_tests(2) = test_case("test_pipeline_reps_not_summing_to_the_replicates", &
                             test_pipeline_reps_not_summing_to_the_replicates)
    all_tests(3) = test_case("test_pipeline_too_few_varying_genes", test_pipeline_too_few_varying_genes)
    all_tests(4) = test_case("test_pipeline_log2_domain", test_pipeline_log2_domain)
    all_tests(5) = test_case("test_pipeline_degree_bounds", test_pipeline_degree_bounds)
    all_tests(6) = test_case("test_pipeline_rejects_nan_and_inf", test_pipeline_rejects_nan_and_inf)
    all_tests(7) = test_case("test_pipeline_dimensions", test_pipeline_dimensions)
  end function get_all_tests_normalization_pipeline

  !> The pipeline is normalize_by_std_dev, optionally quantile_normalization, calc_tiss_avg and
  !| log2_transformation in a row, and must give what those give one by one. Six replicates as one,
  !| two and three tissues: the pipeline borrows output columns 2 and 3 as scratch only where they
  !| exist, so each of those counts takes a different branch.
  subroutine test_pipeline_matches_its_steps()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr, by_std_dev, quantile
    real(real64), allocatable :: from_pipeline(:, :), averages(:, :), from_steps(:, :)
    integer(int32), allocatable :: reps_per_tissue(:)
    real(real64) :: rank_means(n_genes), tmp_row(n_genes)
    integer(int32) :: tmp_perm(n_genes), ierr, n_tissues, i_tissue
    logical(c_bool) :: use_quantile
    character(len=64) :: label

    call fill_linear_trend(expr)
    call normalize_by_std_dev(n_genes, n_replicates, expr, by_std_dev, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_pipeline_matches_its_steps: normalize_by_std_dev")
    call quantile_normalization_expert(n_genes, n_replicates, by_std_dev, quantile, rank_means, tmp_row, tmp_perm, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_pipeline_matches_its_steps: quantile_normalization")

    do n_tissues = 1, 3
      reps_per_tissue = [(n_replicates/n_tissues, i_tissue = 1, n_tissues)]
      allocate (from_pipeline(n_tissues, n_genes), averages(n_tissues, n_genes), from_steps(n_tissues, n_genes))

      do i_tissue = 0, 1
        use_quantile = i_tissue == 1
        write (label, '(a, i0, a, l1)') "test_pipeline_matches_its_steps: tissues=", n_tissues, " quantile=", use_quantile

        if (use_quantile) then
          call calc_tiss_avg(n_genes, n_replicates, n_tissues, reps_per_tissue, quantile, averages, ierr)
        else
          call calc_tiss_avg(n_genes, n_replicates, n_tissues, reps_per_tissue, by_std_dev, averages, ierr)
        end if
        call assert_equal_int(get_err_code(ierr), ERR_OK, trim(label)//" calc_tiss_avg")
        call log2_transformation(n_genes, n_tissues, averages, from_steps, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, trim(label)//" log2_transformation")

        call normalization_pipeline(n_genes, n_replicates, expr, from_pipeline, reps_per_tissue, n_tissues, &
                                    use_quantile=use_quantile, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, trim(label)//" pipeline ierr")
        call assert_equal_array_real(from_pipeline, from_steps, n_tissues*n_genes, 1d-12, trim(label)//" pipeline vs steps")
      end do

      deallocate (from_pipeline, averages, from_steps)
    end do
  end subroutine test_pipeline_matches_its_steps

  !> The replicates per tissue must add up to the rows of `expr`: ERR_INVALID_INPUT otherwise.
  subroutine test_pipeline_reps_not_summing_to_the_replicates()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr
    real(real64), dimension(2, n_genes) :: log_transformed_expr
    integer(int32) :: reps_per_tissue(2), ierr

    call fill_linear_trend(expr)
    reps_per_tissue = [3, 2]
    call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, 2, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "reps [3, 2] sum to 5, but expr has 6 rows")
  end subroutine test_pipeline_reps_not_summing_to_the_replicates

  !> The first step needs five genes that vary; an all-zero matrix has none: ERR_INVALID_INPUT.
  subroutine test_pipeline_too_few_varying_genes()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr
    real(real64), dimension(2, n_genes) :: log_transformed_expr
    integer(int32) :: reps_per_tissue(2), ierr

    expr = 0.0_real64
    reps_per_tissue = [3, 3]
    call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, 2, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_too_few_varying_genes: all zero")
  end subroutine test_pipeline_too_few_varying_genes

  !> The last step takes log2(x + 1), defined only above -1. Negated trend data normalizes to about
  !| -4 .. -7, so the pipeline stops there with ERR_INVALID_INPUT rather than write NaN.
  subroutine test_pipeline_log2_domain()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr
    real(real64), dimension(2, n_genes) :: log_transformed_expr
    integer(int32) :: reps_per_tissue(2), ierr

    call fill_linear_trend(expr)
    expr = -expr
    reps_per_tissue = [3, 3]
    call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, 2, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_log2_domain: averages below -1")
  end subroutine test_pipeline_log2_domain

  !> The pipeline passes span and degree to LOESS, so they follow tox_loess's rules: degree -1 and
  !| 3 are ERR_INVALID_INPUT instead of reaching netlib's ehg182, which stops the whole program, and
  !| so is a span above 1.
  subroutine test_pipeline_degree_bounds()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6, n_tissues = 2
    integer(int32), parameter :: invalid(2) = [-1, 3]
    real(real64), dimension(n_replicates, n_genes) :: expr
    real(real64), dimension(n_tissues, n_genes) :: log_transformed_expr
    integer(int32), dimension(n_tissues) :: reps_per_tissue
    integer(int32) :: ierr, i_degree

    call fill_linear_trend(expr)
    reps_per_tissue = [3, 3]

    do i_degree = 1, size(invalid)
      call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, &
                                  degree=invalid(i_degree), ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_degree_bounds: degrees -1 and 3 are invalid")
    end do

    call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, n_tissues, &
                                span=1.5_real64, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_degree_bounds: a span above 1 is invalid")
  end subroutine test_pipeline_degree_bounds

  !> NaN and Inf are rejected up front, before the first step: TOX writes no NaN.
  subroutine test_pipeline_rejects_nan_and_inf()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr
    real(real64), dimension(2, n_genes) :: log_transformed_expr
    real(real64) :: bad(2)
    integer(int32) :: reps_per_tissue(2), ierr, i_bad

    call fill_linear_trend(expr)
    reps_per_tissue = [3, 3]
    bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    do i_bad = 1, size(bad)
      expr(1, 1) = bad(i_bad)
      call normalization_pipeline(n_genes, n_replicates, expr, log_transformed_expr, reps_per_tissue, 2, ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                            "test_pipeline_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
    end do
  end subroutine test_pipeline_rejects_nan_and_inf

  !> Each dimension on its own. n_genes = 0 is ERR_EMPTY_INPUT and a negative one ERR_INVALID_INPUT.
  !| n_replicates is checked against sum(reps_per_tissue) instead of as a dimension (DM_MIN/DM_MAX,
  !| until #203), so n_replicates = 0 is ERR_INVALID_INPUT, and so is n_tissues = 0: its empty sum
  !| is 0, which n_replicates misses before n_tissues' own check runs.
  subroutine test_pipeline_dimensions()
    real(real64), dimension(1, 1) :: expr, log_transformed_expr
    integer(int32) :: reps_per_tissue(1), ierr

    expr = 1.0_real64
    reps_per_tissue = [1]
    call normalization_pipeline(0, 1, expr, log_transformed_expr, reps_per_tissue, 1, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_pipeline_dimensions: n_genes = 0")
    call normalization_pipeline(1, 0, expr, log_transformed_expr, reps_per_tissue, 1, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_dimensions: n_replicates = 0")
    call normalization_pipeline(1, 1, expr, log_transformed_expr, reps_per_tissue, 0, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_dimensions: n_tissues = 0")
    call normalization_pipeline(-1, 1, expr, log_transformed_expr, reps_per_tissue, 1, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_pipeline_dimensions: n_genes = -1")
  end subroutine test_pipeline_dimensions

end module mod_test_normalization_pipeline
