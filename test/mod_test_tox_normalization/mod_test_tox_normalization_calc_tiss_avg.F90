!> The `calc_tiss_avg` cases: hand-derived tissue averages for every way replicates can be grouped,
!| the bounds on the replicate counts, and the input checks.
module mod_test_tox_normalization_calc_tiss_avg
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use tox_normalization
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_tox_normalization_calc_tiss_avg() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(10))

    all_tests(1) = test_case("test_calc_tiss_avg_three_tissues", test_calc_tiss_avg_three_tissues)
    all_tests(2) = test_case("test_calc_tiss_avg_single_tissue", test_calc_tiss_avg_single_tissue)
    all_tests(3) = test_case("test_calc_tiss_avg_unequal_replicates", test_calc_tiss_avg_unequal_replicates)
    all_tests(4) = test_case("test_calc_tiss_avg_single_replicate", test_calc_tiss_avg_single_replicate)
    all_tests(5) = test_case("test_calc_tiss_avg_magnitudes_and_signs", test_calc_tiss_avg_magnitudes_and_signs)
    all_tests(6) = test_case("test_calc_tiss_avg_reps_bounds", test_calc_tiss_avg_reps_bounds)
    all_tests(7) = test_case("test_calc_tiss_avg_reps_not_summing_to_the_replicates", &
                             test_calc_tiss_avg_reps_not_summing_to_the_replicates)
    all_tests(8) = test_case("test_calc_tiss_avg_dimensions", test_calc_tiss_avg_dimensions)
    all_tests(9) = test_case("test_calc_tiss_avg_rejects_nan_and_inf", test_calc_tiss_avg_rejects_nan_and_inf)
    all_tests(10) = test_case("test_calc_tiss_avg_extreme_magnitudes", test_calc_tiss_avg_extreme_magnitudes)
  end function get_all_tests_tox_normalization_calc_tiss_avg

  !> An average is never larger than its largest value, so it must not overflow where the sum
  !| does: two replicates of huge average to huge, and [huge, huge/2] to 0.75*huge. Summing first
  !| overflowed both to Inf.
  subroutine test_calc_tiss_avg_extreme_magnitudes()
    integer(int32) :: ierr
    integer(int32), dimension(1) :: reps_per_tissue
    real(real64), dimension(2, 2) :: expr
    real(real64), dimension(1, 2) :: averages, expected

    expr(:, 1) = [huge(1d0), huge(1d0)]
    expr(:, 2) = [huge(1d0), 0.5d0*huge(1d0)]
    reps_per_tissue = [2]
    expected(1, :) = [huge(1d0), 0.75d0*huge(1d0)]

    call calc_tiss_avg(2, 2, 1, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_extreme_magnitudes: ierr")
    call assert_equal_array_real(averages, expected, 2, 0d0, "test_calc_tiss_avg_extreme_magnitudes: averages of huge values")
  end subroutine test_calc_tiss_avg_extreme_magnitudes

  !> Three tissues of two replicates: rows 1-2, 3-4 and 5-6 of each gene are averaged.
  subroutine test_calc_tiss_avg_three_tissues()
    integer(int32) :: ierr
    integer(int32), dimension(3) :: reps_per_tissue
    real(real64), dimension(6, 2) :: expr
    real(real64), dimension(3, 2) :: averages, expected

    expr(:, 1) = [1d0, 3d0, 5d0, 2d0, 4d0, 6d0]
    expr(:, 2) = [7d0, 9d0, 11d0, 8d0, 10d0, 12d0]
    reps_per_tissue = [2, 2, 2]
    ! gene 1: mean(1, 3) = 2, mean(5, 2) = 3.5, mean(4, 6) = 5; gene 2: 8, 9.5, 11
    expected(:, 1) = [2d0, 3.5d0, 5d0]
    expected(:, 2) = [8d0, 9.5d0, 11d0]

    call calc_tiss_avg(2, 6, 3, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_three_tissues: ierr")
    call assert_equal_array_real(averages, expected, 6, 0d0, "test_calc_tiss_avg_three_tissues: averages")
  end subroutine test_calc_tiss_avg_three_tissues

  !> One tissue holding every replicate averages each whole gene.
  subroutine test_calc_tiss_avg_single_tissue()
    integer(int32) :: ierr
    integer(int32), dimension(1) :: reps_per_tissue
    real(real64), dimension(2, 3) :: expr
    real(real64), dimension(1, 3) :: averages, expected

    expr(:, 1) = [1d0, 4d0]
    expr(:, 2) = [2d0, 5d0]
    expr(:, 3) = [3d0, 6d0]
    reps_per_tissue = [2]
    expected(1, :) = [2.5d0, 3.5d0, 4.5d0]

    call calc_tiss_avg(3, 2, 1, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_single_tissue: ierr")
    call assert_equal_array_real(averages, expected, 3, 0d0, "test_calc_tiss_avg_single_tissue: averages")
  end subroutine test_calc_tiss_avg_single_tissue

  !> Tissues of 2, 3 and 2 replicates: the group boundaries move with the counts.
  subroutine test_calc_tiss_avg_unequal_replicates()
    integer(int32) :: ierr
    integer(int32), dimension(3) :: reps_per_tissue
    real(real64), dimension(7, 2) :: expr
    real(real64), dimension(3, 2) :: averages, expected

    expr(:, 1) = [1d0, 3d0, 5d0, 7d0, 9d0, 11d0, 13d0]
    expr(:, 2) = [2d0, 4d0, 6d0, 8d0, 10d0, 12d0, 14d0]
    reps_per_tissue = [2, 3, 2]
    ! gene 1: mean(1, 3) = 2, mean(5, 7, 9) = 7, mean(11, 13) = 12; gene 2: 3, 8, 13
    expected(:, 1) = [2d0, 7d0, 12d0]
    expected(:, 2) = [3d0, 8d0, 13d0]

    call calc_tiss_avg(2, 7, 3, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_unequal_replicates: ierr")
    ! A few ulps: an optimizer may divide by 3 as a multiplication with the rounded 1/3 (ifx does),
    ! which the tissues of 2 replicates never see, 1/2 being exact.
    call assert_equal_array_real(averages, expected, 6, 1d-14, "test_calc_tiss_avg_unequal_replicates: averages")
  end subroutine test_calc_tiss_avg_unequal_replicates

  !> One replicate per tissue: the averages are the input itself.
  subroutine test_calc_tiss_avg_single_replicate()
    integer(int32) :: ierr
    integer(int32), dimension(3) :: reps_per_tissue
    real(real64), dimension(3, 2) :: expr, averages

    expr(:, 1) = [1d0, 2d0, 3d0]
    expr(:, 2) = [4d0, 5d0, 6d0]
    reps_per_tissue = [1, 1, 1]

    call calc_tiss_avg(2, 3, 3, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_single_replicate: ierr")
    call assert_equal_array_real(averages, expr, 6, 0d0, "test_calc_tiss_avg_single_replicate: must equal the input")
  end subroutine test_calc_tiss_avg_single_replicate

  !> Averages are exact across fifteen orders of magnitude, and signs cancel: all of these are
  !| exactly representable, so no tolerance is needed.
  subroutine test_calc_tiss_avg_magnitudes_and_signs()
    integer(int32) :: ierr
    integer(int32), dimension(2) :: reps_per_tissue
    real(real64), dimension(4, 3) :: expr
    real(real64), dimension(2, 3) :: averages, expected

    expr(:, 1) = [1d6, 1d9, 1d12, 1d15]
    expr(:, 2) = [2d3, 2d9, 2d9, 2d15]
    expr(:, 3) = [-1d0, -3d0, -7d0, 7d0]
    reps_per_tissue = [2, 2]
    expected(:, 1) = [5.005d8, 5.005d14]
    expected(:, 2) = [1.000001d9, 1.000001d15]
    expected(:, 3) = [-2d0, 0d0]

    call calc_tiss_avg(3, 4, 2, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_tiss_avg_magnitudes_and_signs: ierr")
    call assert_equal_array_real(averages, expected, 6, 0d0, "test_calc_tiss_avg_magnitudes_and_signs: averages")
  end subroutine test_calc_tiss_avg_magnitudes_and_signs

  !> Every tissue needs at least one replicate: a count of 0 or -1 is ERR_INVALID_INPUT.
  subroutine test_calc_tiss_avg_reps_bounds()
    integer(int32) :: ierr
    integer(int32), dimension(2) :: reps_per_tissue
    real(real64), dimension(2, 1) :: expr
    real(real64), dimension(2, 1) :: averages

    expr = 1d0

    reps_per_tissue = [2, 0]
    call calc_tiss_avg(1, 2, 2, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_tiss_avg_reps_bounds: a tissue with 0 replicates")

    reps_per_tissue = [-1, 3]
    call calc_tiss_avg(1, 2, 2, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_tiss_avg_reps_bounds: a tissue with -1 replicates")
  end subroutine test_calc_tiss_avg_reps_bounds

  !> The replicates must add up to the rows of `expr`. Fewer would shift every gene after the
  !| first onto the wrong rows; more would read past the end of the matrix.
  subroutine test_calc_tiss_avg_reps_not_summing_to_the_replicates()
    integer(int32) :: ierr
    integer(int32), dimension(2) :: reps_per_tissue
    real(real64), dimension(6, 2) :: expr
    real(real64), dimension(2, 2) :: averages

    expr = 1d0

    reps_per_tissue = [3, 2]
    call calc_tiss_avg(2, 6, 2, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "reps [3, 2] sum to 5, but expr has 6 rows")
    reps_per_tissue = [3, 4]
    call calc_tiss_avg(2, 6, 2, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "reps [3, 4] sum to 7, but expr has 6 rows")
  end subroutine test_calc_tiss_avg_reps_not_summing_to_the_replicates

  !> Each dimension on its own. n_genes = 0 is ERR_EMPTY_INPUT and a negative one ERR_INVALID_INPUT.
  !| n_replicates is checked against sum(reps_per_tissue) instead of as a dimension (DM_MIN/DM_MAX,
  !| until #203), so n_replicates = 0 is ERR_INVALID_INPUT, and so is n_tissues = 0: its empty sum
  !| is 0, which n_replicates misses before n_tissues' own check runs.
  subroutine test_calc_tiss_avg_dimensions()
    integer(int32) :: ierr
    integer(int32), dimension(1) :: reps_per_tissue
    real(real64), dimension(1, 1) :: expr, averages

    expr = 1d0
    reps_per_tissue = [1]

    call calc_tiss_avg(0, 1, 1, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_calc_tiss_avg_dimensions: n_genes = 0")
    call calc_tiss_avg(1, 0, 1, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_tiss_avg_dimensions: n_replicates = 0")
    call calc_tiss_avg(1, 1, 0, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_tiss_avg_dimensions: n_tissues = 0")
    call calc_tiss_avg(-1, 1, 1, reps_per_tissue, expr, averages, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_tiss_avg_dimensions: n_genes = -1")
  end subroutine test_calc_tiss_avg_dimensions

  !> NaN and Inf are rejected up front rather than carried into the result: TOX writes no NaN.
  subroutine test_calc_tiss_avg_rejects_nan_and_inf()
    integer(int32) :: ierr, i_bad
    integer(int32), dimension(2) :: reps_per_tissue
    real(real64), dimension(6, 2) :: expr
    real(real64), dimension(2, 2) :: averages
    real(real64) :: bad(2)

    expr = 1d0
    reps_per_tissue = [3, 3]
    bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    do i_bad = 1, size(bad)
      expr(1, 1) = bad(i_bad)
      call calc_tiss_avg(2, 6, 2, reps_per_tissue, expr, averages, ierr)
      call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                            "test_calc_tiss_avg_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
    end do
  end subroutine test_calc_tiss_avg_rejects_nan_and_inf

end module mod_test_tox_normalization_calc_tiss_avg
