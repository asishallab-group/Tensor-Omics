!> The `normalize_by_std_dev` cases. Exact values need data LOESS reproduces exactly: genes on a
!| linear mean-sd trend (mod_test_tox_normalization_fixtures), whose normalization has a closed
!| form. Around that: an off-trend gene, genes without variance, the bounds, and the input checks.
module mod_test_tox_normalization_normalize_by_std_dev
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use tox_normalization
  ! the tox_normalization module used to re-export it; it is f42 infrastructure
  use f42_math_impl, only: std_dev
  use mod_test_tox_normalization_fixtures, only: fill_linear_trend, linear_trend_offset, linear_trend_sd
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public

  ! LOESS solves a least-squares system: a straight line comes back exact to rounding, not bits
  real(real64), parameter :: TOL = 1.0e-9_real64

contains

  !> Get array of all available tests.
  function get_all_tests_tox_normalization_normalize_by_std_dev() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(10))

    all_tests(1) = test_case("test_std_dev", test_std_dev)
    all_tests(2) = test_case("test_std_dev_linear_trend", test_std_dev_linear_trend)
    all_tests(3) = test_case("test_std_dev_off_trend_gene", test_std_dev_off_trend_gene)
    all_tests(4) = test_case("test_std_dev_zero_variance_genes", test_std_dev_zero_variance_genes)
    all_tests(5) = test_case("test_std_dev_needs_five_varying_genes", test_std_dev_needs_five_varying_genes)
    all_tests(6) = test_case("test_std_dev_span_bounds", test_std_dev_span_bounds)
    all_tests(7) = test_case("test_std_dev_degree_bounds", test_std_dev_degree_bounds)
    all_tests(8) = test_case("test_std_dev_rejects_nan_and_inf", test_std_dev_rejects_nan_and_inf)
    all_tests(9) = test_case("test_std_dev_dimensions", test_std_dev_dimensions)
    all_tests(10) = test_case("test_std_dev_negative_fit_keeps_sign", test_std_dev_negative_fit_keeps_sign)
  end function get_all_tests_tox_normalization_normalize_by_std_dev

  !> A LOESS fit can dip below zero on non-negative data, and dividing a gene by a negative fitted
  !| sd flips its sign. Fifteen genes whose sd stays near 0 and then climbs steeply, fitted with the
  !| default span and degree, do exactly that for genes 4 to 6. Every value here is positive, so
  !| every normalized value must be too: a non-positive fit falls back to the gene's own sd, as a
  !| fit near zero already did.
  subroutine test_std_dev_negative_fit_keeps_sign()
    integer(int32), parameter :: n_genes = 15, n_replicates = 2
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized
    real(real64) :: gene_sd(n_genes), gene_mean
    integer(int32) :: ierr, i_gene

    ! gene i: mean 100 + 10*i, replicates mean -/+ sd, so a population sd of gene_sd(i)
    gene_sd = [1d-3, 1d-3, 1d-3, 1d-3, 1d-3, 1d-3, 1d-3, 0.01d0, 0.05d0, 0.2d0, 0.8d0, 2.5d0, 6d0, 12d0, 20d0]
    do i_gene = 1, n_genes
      gene_mean = 100d0 + 10d0*real(i_gene, real64)
      expr(:, i_gene) = [gene_mean - gene_sd(i_gene), gene_mean + gene_sd(i_gene)]
    end do

    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_negative_fit_keeps_sign: ierr")
    call assert_true(all(normalized > 0d0), "test_std_dev_negative_fit_keeps_sign: a positive gene must stay positive")
  end subroutine test_std_dev_negative_fit_keeps_sign

  !> f42_math's std_dev, which the procedure builds its trend from. It belongs to the f42 suites
  !| and moves there with the f42 batch.
  subroutine test_std_dev()
    real(real64) :: v(5), w(1)

    ! population variance of 1..5 is 2, sample variance 2.5
    v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    call assert_equal_real(std_dev(v), sqrt(2.0_real64), 1e-12_real64, "test_std_dev: population std_dev")
    call assert_equal_real(std_dev(v, do_bessel_correction=.true._c_bool), sqrt(2.5_real64), 1e-12_real64, &
                           "test_std_dev: sample std_dev (Bessel)")

    ! a single value has no spread in either mode
    w = [3.14159_real64]
    call assert_equal_real(std_dev(w), 0.0_real64, 0.0_real64, "test_std_dev: length 1 (population)")
    call assert_equal_real(std_dev(w, do_bessel_correction=.true._c_bool), 0.0_real64, 0.0_real64, &
                           "test_std_dev: length 1 (sample)")
  end subroutine test_std_dev

  !> Genes exactly on the linear trend are divided by their own sd, as LOESS reproduces the line:
  !| gene i, i*(10 + d_j), becomes (10 + d_j)/s for every gene.
  subroutine test_std_dev_linear_trend()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
    integer(int32) :: ierr, i_replicate

    call fill_linear_trend(expr)
    do i_replicate = 1, n_replicates
      expected(i_replicate, :) = (10.0_real64 + linear_trend_offset(i_replicate, n_replicates))/linear_trend_sd(n_replicates)
    end do

    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_linear_trend: ierr")
    call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                 "test_std_dev_linear_trend: every gene becomes (10 + d_j)/s")
  end subroutine test_std_dev_linear_trend

  !> The point of fitting a trend: a gene off it is divided by the trend's sd at its mean, not by
  !| its own. Gene 10 keeps its mean of 100 but has three times the trend's spread, 10*(10 + 3*d_j);
  !| the robust fit discounts it, so it becomes (10 + 3*d_j)/s while the others stay (10 + d_j)/s.
  subroutine test_std_dev_off_trend_gene()
    integer(int32), parameter :: n_genes = 20, n_replicates = 6, off_trend = 10
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
    real(real64) :: offset, s
    integer(int32) :: ierr, i_replicate

    call fill_linear_trend(expr)
    s = linear_trend_sd(n_replicates)
    do i_replicate = 1, n_replicates
      offset = linear_trend_offset(i_replicate, n_replicates)
      expr(i_replicate, off_trend) = real(off_trend, real64)*(10.0_real64 + 3.0_real64*offset)
      expected(i_replicate, :) = (10.0_real64 + offset)/s
      expected(i_replicate, off_trend) = (10.0_real64 + 3.0_real64*offset)/s
    end do

    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_off_trend_gene: ierr")
    call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                 "test_std_dev_off_trend_gene: gene 10 must be divided by the trend's sd, not its own")
  end subroutine test_std_dev_off_trend_gene

  !> Genes without variance say nothing about the trend: they are left out of the fit and returned
  !| unchanged, while the genes on the trend are normalized as usual.
  subroutine test_std_dev_zero_variance_genes()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6, n_varying = 7
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized, expected
    integer(int32) :: ierr, i_replicate

    call fill_linear_trend(expr(:, 1:n_varying))
    expr(:, n_varying + 1:) = 1.0_real64
    expected(:, n_varying + 1:) = 1.0_real64
    do i_replicate = 1, n_replicates
      expected(i_replicate, 1:n_varying) = (10.0_real64 + linear_trend_offset(i_replicate, n_replicates)) &
                                           /linear_trend_sd(n_replicates)
    end do

    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_zero_variance_genes: ierr")
    call assert_equal_array_real(normalized, expected, n_genes*n_replicates, TOL, &
                                 "test_std_dev_zero_variance_genes: constant genes unchanged, the rest normalized")
  end subroutine test_std_dev_zero_variance_genes

  !> The fit needs five genes that vary: four is ERR_INVALID_INPUT, five is enough. The span is 1
  !| here, so that five points are also enough for LOESS itself; the bound under test is the count.
  subroutine test_std_dev_needs_five_varying_genes()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized
    integer(int32) :: ierr

    call fill_linear_trend(expr(:, 1:4))
    expr(:, 5:) = 1.0_real64
    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, span=1.0_real64, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_std_dev_needs_five_varying_genes: four varying genes")

    call fill_linear_trend(expr(:, 1:5))
    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, span=1.0_real64, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_needs_five_varying_genes: five varying genes")
  end subroutine test_std_dev_needs_five_varying_genes

  !> The span follows tox_loess's rules, as it is the LOESS span: at most 1, so 1 is valid and 1.5
  !| and 3 are ERR_INVALID_INPUT. At the low end 0, a negative span, and 0.1 of ten genes (too few
  !| points for LOESS) are ERR_INVALID_INPUT. A NaN span is ERR_NAN_INF.
  subroutine test_std_dev_span_bounds()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64), parameter :: invalid(5) = [0.0_real64, -1.0_real64, 0.1_real64, 1.5_real64, 3.0_real64]
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized
    integer(int32) :: ierr, i_span

    call fill_linear_trend(expr)
    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, span=1.0_real64, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_span_bounds: span 1 is valid")
    do i_span = 1, size(invalid)
      call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, span=invalid(i_span), ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_std_dev_span_bounds: span outside (0, 1]")
    end do

    call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, span=ieee_value(1.0_real64, ieee_quiet_nan), ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, "test_std_dev_span_bounds: NaN span")
  end subroutine test_std_dev_span_bounds

  !> degree must be 0, 1 or 2, the local models LOESS knows; anything else is ERR_INVALID_INPUT.
  !| Unchecked, -1 and 3 reached netlib's ehg182, which stops the whole program -- the caller's
  !| Python or R session with it.
  subroutine test_std_dev_degree_bounds()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    integer(int32), parameter :: valid(3) = [0, 1, 2], invalid(2) = [-1, 3]
    real(real64), dimension(n_replicates, n_genes) :: expr, normalized
    integer(int32) :: ierr, i_degree

    call fill_linear_trend(expr)
    do i_degree = 1, size(valid)
      call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, degree=valid(i_degree), ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_OK, "test_std_dev_degree_bounds: degrees 0, 1 and 2 are valid")
    end do
    do i_degree = 1, size(invalid)
      call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, degree=invalid(i_degree), ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_std_dev_degree_bounds: degrees -1 and 3 are invalid")
    end do
  end subroutine test_std_dev_degree_bounds

  !> NaN and Inf are rejected up front rather than carried into the result: TOX writes no NaN.
  !| The matrix varies, so without that check the call would succeed and hide the NaN.
  subroutine test_std_dev_rejects_nan_and_inf()
    integer(int32), parameter :: n_genes = 10, n_replicates = 6
    real(real64) :: expr(n_replicates, n_genes), normalized(n_replicates, n_genes), bad(2)
    integer(int32) :: ierr, i_bad

    call fill_linear_trend(expr)
    bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    do i_bad = 1, size(bad)
      expr(1, 1) = bad(i_bad)
      call normalize_by_std_dev(n_genes, n_replicates, expr, normalized, ierr=ierr)
      call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                            "test_std_dev_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
    end do
  end subroutine test_std_dev_rejects_nan_and_inf

  !> Each dimension on its own: zero is ERR_EMPTY_INPUT, negative ERR_INVALID_INPUT.
  subroutine test_std_dev_dimensions()
    real(real64), dimension(1, 1) :: expr, normalized
    integer(int32) :: ierr

    expr = 1.0_real64
    call normalize_by_std_dev(0, 1, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_std_dev_dimensions: n_genes = 0")
    call normalize_by_std_dev(1, 0, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_std_dev_dimensions: n_replicates = 0")
    call normalize_by_std_dev(-1, 1, expr, normalized, ierr=ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_std_dev_dimensions: n_genes = -1")
  end subroutine test_std_dev_dimensions

end module mod_test_tox_normalization_normalize_by_std_dev
