!> The `calc_fchange` cases: hand-derived fold changes for every way a pair can be formed, the
!| index bounds, and the input checks.
module mod_test_tox_normalization_calc_fchange
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
  function get_all_tests_tox_normalization_calc_fchange() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(4))

    all_tests(1) = test_case("test_calc_fchange_values", test_calc_fchange_values)
    all_tests(2) = test_case("test_calc_fchange_index_bounds", test_calc_fchange_index_bounds)
    all_tests(3) = test_case("test_calc_fchange_dimensions", test_calc_fchange_dimensions)
    all_tests(4) = test_case("test_calc_fchange_rejects_nan_and_inf", test_calc_fchange_rejects_nan_and_inf)
  end function get_all_tests_tox_normalization_calc_fchange

  !> Each pair is condition minus control, gene by gene. Three tissues, two genes, and four pairs
  !| covering what a caller can do: two conditions against one control, a condition listed before
  !| its control, and a tissue against itself (which must give exactly 0).
  subroutine test_calc_fchange_values()
    integer(int32), parameter :: n_genes = 2, n_tissues = 3, n_pairs = 4
    integer(int32), dimension(n_pairs) :: control_tissues, condition_tissues
    real(real64), dimension(n_tissues, n_genes) :: expr
    real(real64), dimension(n_pairs, n_genes) :: fold_changes, expected
    integer(int32) :: ierr

    ! tissue rows: tissue 1 = [1, 10], tissue 2 = [4, 7], tissue 3 = [9, 15]
    expr(:, 1) = [1.0_real64, 4.0_real64, 9.0_real64]
    expr(:, 2) = [10.0_real64, 7.0_real64, 15.0_real64]
    control_tissues = [1, 1, 3, 2]
    condition_tissues = [2, 3, 2, 2]
    expected(1, :) = [3.0_real64, -3.0_real64]   ! 2 - 1: [4 - 1, 7 - 10]
    expected(2, :) = [8.0_real64, 5.0_real64]    ! 3 - 1: [9 - 1, 15 - 10]
    expected(3, :) = [-5.0_real64, -8.0_real64]  ! 2 - 3: [4 - 9, 7 - 15]
    expected(4, :) = [0.0_real64, 0.0_real64]    ! 2 - 2

    call calc_fchange(n_genes, n_tissues, n_pairs, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_fchange_values: ierr")
    call assert_equal_array_real(fold_changes, expected, n_pairs*n_genes, 0.0_real64, &
                                 "test_calc_fchange_values: condition minus control")
  end subroutine test_calc_fchange_values

  !> A tissue index must lie in 1..n_tissues, for the control and the condition alike. The last
  !| tissue is valid; 0, -1 and n_tissues + 1 are ERR_INVALID_INPUT.
  subroutine test_calc_fchange_index_bounds()
    integer(int32), parameter :: n_genes = 2, n_tissues = 3
    integer(int32), dimension(1) :: control_tissues, condition_tissues
    real(real64), dimension(n_tissues, n_genes) :: expr
    real(real64), dimension(1, n_genes) :: fold_changes
    integer(int32) :: ierr, i_bad
    integer(int32), parameter :: out_of_range(3) = [0, -1, n_tissues + 1]

    expr = 1.0_real64

    control_tissues = [1]
    condition_tissues = [n_tissues]
    call calc_fchange(n_genes, n_tissues, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_fchange_index_bounds: the last tissue is valid")

    do i_bad = 1, size(out_of_range)
      control_tissues = [out_of_range(i_bad)]
      condition_tissues = [1]
      call calc_fchange(n_genes, n_tissues, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
      call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_fchange_index_bounds: control out of range")

      control_tissues = [1]
      condition_tissues = [out_of_range(i_bad)]
      call calc_fchange(n_genes, n_tissues, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
      call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_fchange_index_bounds: condition out of range")
    end do
  end subroutine test_calc_fchange_index_bounds

  !> Each dimension on its own: zero is ERR_EMPTY_INPUT, negative ERR_INVALID_INPUT.
  subroutine test_calc_fchange_dimensions()
    integer(int32), dimension(1) :: control_tissues, condition_tissues
    real(real64), dimension(1, 1) :: expr, fold_changes
    integer(int32) :: ierr

    expr = 1.0_real64
    control_tissues = [1]
    condition_tissues = [1]

    call calc_fchange(0, 1, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_calc_fchange_dimensions: n_genes = 0")
    call calc_fchange(1, 0, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_calc_fchange_dimensions: n_tissues = 0")
    call calc_fchange(1, 1, 0, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_calc_fchange_dimensions: n_pairs = 0")
    call calc_fchange(-1, 1, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
    call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_calc_fchange_dimensions: n_genes = -1")
  end subroutine test_calc_fchange_dimensions

  !> NaN and Inf are rejected up front rather than carried into the result: TOX writes no NaN.
  subroutine test_calc_fchange_rejects_nan_and_inf()
    integer(int32) :: ierr, i_bad
    integer(int32), dimension(1) :: control_tissues, condition_tissues
    real(real64), dimension(3, 2) :: expr
    real(real64), dimension(1, 2) :: fold_changes
    real(real64) :: bad(2)

    expr = 1d0
    control_tissues = [1]
    condition_tissues = [2]
    bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    do i_bad = 1, size(bad)
      expr(2, 1) = bad(i_bad)
      call calc_fchange(2, 3, 1, control_tissues, condition_tissues, expr, fold_changes, ierr)
      call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                            "test_calc_fchange_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
    end do
  end subroutine test_calc_fchange_rejects_nan_and_inf

end module mod_test_tox_normalization_calc_fchange
