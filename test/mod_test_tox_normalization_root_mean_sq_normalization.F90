!> The `root_mean_sq_normalization` cases: hand-derived RMS scalings, the properties the scaling
!| has, the genes it leaves alone, magnitudes past the real64 range, and the input checks.
module mod_test_tox_normalization_root_mean_sq_normalization
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use tox_normalization
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public

  real(real64), parameter :: TOL = 4*epsilon(1.0_real64)

contains

  !> Get array of all available tests.
  function get_all_tests_tox_normalization_root_mean_sq_normalization() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(8))

    all_tests(1) = test_case("test_rms_values", test_rms_values)
    all_tests(2) = test_case("test_rms_identity", test_rms_identity)
    all_tests(3) = test_case("test_rms_scale_and_sign", test_rms_scale_and_sign)
    all_tests(4) = test_case("test_rms_near_zero_genes_pass_through", test_rms_near_zero_genes_pass_through)
    all_tests(5) = test_case("test_rms_single_replicate", test_rms_single_replicate)
    all_tests(6) = test_case("test_rms_extreme_magnitudes", test_rms_extreme_magnitudes)
    all_tests(7) = test_case("test_rms_rejects_nan_and_inf", test_rms_rejects_nan_and_inf)
    all_tests(8) = test_case("test_rms_dimensions", test_rms_dimensions)
  end function get_all_tests_tox_normalization_root_mean_sq_normalization

  !> Each gene is divided by sqrt(mean(x**2)) over its replicates: [1, 7] has RMS 5 and becomes
  !| [0.2, 1.4], [-3, 3] has RMS 3 and becomes [-1, 1], and a constant [5, 5] becomes [1, 1].
  subroutine test_rms_values()
    real(real64), dimension(2, 3) :: expr, normalized, expected
    integer(int32) :: ierr

    expr(:, 1) = [1.0_real64, 7.0_real64]
    expr(:, 2) = [-3.0_real64, 3.0_real64]
    expr(:, 3) = [5.0_real64, 5.0_real64]
    expected(:, 1) = [0.2_real64, 1.4_real64]
    expected(:, 2) = [-1.0_real64, 1.0_real64]
    expected(:, 3) = [1.0_real64, 1.0_real64]

    call root_mean_sq_normalization(3, 2, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_values: ierr")
    call assert_equal_array_real(normalized, expected, 6, TOL, "test_rms_values: x/sqrt(mean(x**2))")
  end subroutine test_rms_values

  !> The identity matrix: each gene has one 1 among three replicates, RMS sqrt(1/3), so its 1
  !| becomes sqrt(3) and its zeros stay exactly 0.
  subroutine test_rms_identity()
    real(real64), dimension(3, 3) :: expr, normalized, expected
    integer(int32) :: ierr, i_gene

    expr = 0.0_real64
    expected = 0.0_real64
    do i_gene = 1, 3
      expr(i_gene, i_gene) = 1.0_real64
      expected(i_gene, i_gene) = sqrt(3.0_real64)
    end do

    call root_mean_sq_normalization(3, 3, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_identity: ierr")
    call assert_equal_array_real(normalized, expected, 9, TOL, "test_rms_identity: diagonal sqrt(3), zeros stay 0")
  end subroutine test_rms_identity

  !> The scaling ignores magnitude and keeps sign: a million times a gene normalizes to the gene's
  !| own result, and its negative to the negated result.
  subroutine test_rms_scale_and_sign()
    real(real64), dimension(3, 3) :: expr, normalized
    real(real64), dimension(3) :: reference
    integer(int32) :: ierr

    expr(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
    expr(:, 2) = 1.0e6_real64*expr(:, 1)
    expr(:, 3) = -expr(:, 1)

    call root_mean_sq_normalization(3, 3, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_scale_and_sign: ierr")
    reference = normalized(:, 1)
    call assert_equal_array_real(normalized(:, 2), reference, 3, TOL, "test_rms_scale_and_sign: 1e6 times the gene")
    reference = -reference
    call assert_equal_array_real(normalized(:, 3), reference, 3, TOL, "test_rms_scale_and_sign: the negated gene")
  end subroutine test_rms_scale_and_sign

  !> A gene whose RMS is at most 1e-12 is left as it is rather than divided by next to nothing: an
  !| all-zero gene stays zero, and so does a gene of 1e-13 values (FES decided to keep this).
  subroutine test_rms_near_zero_genes_pass_through()
    real(real64), dimension(2, 2) :: expr, normalized
    integer(int32) :: ierr

    expr(:, 1) = 0.0_real64
    expr(:, 2) = [1.0e-13_real64, -1.0e-13_real64]

    call root_mean_sq_normalization(2, 2, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_near_zero_genes_pass_through: ierr")
    call assert_equal_array_real(normalized, expr, 4, 0.0_real64, "test_rms_near_zero_genes_pass_through: must stay unchanged")
  end subroutine test_rms_near_zero_genes_pass_through

  !> With one replicate a gene's RMS is its own magnitude, so every gene becomes its sign.
  subroutine test_rms_single_replicate()
    real(real64), dimension(1, 2) :: expr, normalized, expected
    integer(int32) :: ierr

    expr(1, :) = [-5.0_real64, 2.0_real64]
    expected(1, :) = [-1.0_real64, 1.0_real64]

    call root_mean_sq_normalization(2, 1, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_single_replicate: ierr")
    call assert_equal_array_real(normalized, expected, 2, 0.0_real64, "test_rms_single_replicate: each gene becomes its sign")
  end subroutine test_rms_single_replicate

  !> A gene whose squares leave the real64 range is still scaled by its RMS: [1e200, 7e200] has
  !| RMS sqrt((1 + 49)/2) * 1e200 = 5e200 and becomes [0.2, 1.4]. Summing the squares directly
  !| overflows to Inf, and dividing by Inf zeroed the gene.
  subroutine test_rms_extreme_magnitudes()
    real(real64) :: expr(2, 1), normalized(2, 1), expected(2, 1)
    integer(int32) :: ierr

    expr(:, 1) = [1.0e200_real64, 7.0e200_real64]
    expected(:, 1) = [0.2_real64, 1.4_real64]

    call root_mean_sq_normalization(1, 2, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_rms_extreme_magnitudes: ierr")
    call assert_equal_array_real(normalized, expected, 2, TOL, &
                                 "test_rms_extreme_magnitudes: [1e200, 7e200] must become [0.2, 1.4]")
  end subroutine test_rms_extreme_magnitudes

  !> NaN and Inf are rejected up front rather than carried into the result: TOX writes no NaN.
  subroutine test_rms_rejects_nan_and_inf()
    integer(int32), parameter :: n_genes = 3, n_replicates = 2
    real(real64) :: expr(n_replicates, n_genes), normalized(n_replicates, n_genes), bad(2)
    integer(int32) :: ierr, i_bad

    expr = 1.0_real64
    bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    do i_bad = 1, size(bad)
      expr(1, 1) = bad(i_bad)
      call root_mean_sq_normalization(n_genes, n_replicates, expr, normalized, ierr)
      call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                            "test_rms_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
    end do
  end subroutine test_rms_rejects_nan_and_inf

  !> Each dimension on its own: zero is ERR_EMPTY_INPUT, negative ERR_INVALID_INPUT.
  subroutine test_rms_dimensions()
    real(real64), dimension(1, 1) :: expr, normalized
    integer(int32) :: ierr

    expr = 1.0_real64
    call root_mean_sq_normalization(0, 1, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_rms_dimensions: n_genes = 0")
    call root_mean_sq_normalization(1, 0, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_rms_dimensions: n_replicates = 0")
    call root_mean_sq_normalization(-1, 1, expr, normalized, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_rms_dimensions: n_genes = -1")
  end subroutine test_rms_dimensions

end module mod_test_tox_normalization_root_mean_sq_normalization
