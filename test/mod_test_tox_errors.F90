!> Unit test suite for the array validators in tox_errors.
!> They answer from a fast scan when an array is clean and from a per-element loop otherwise, so every
!> test here puts the problem at every position of an array long enough to span several vector
!> lengths: the vector body, its remainder and the seam between them must all see it. This is also
!> the guard against a compiler assuming NaN cannot occur, so it matters under every build profile.

module mod_test_tox_errors
  use asserts
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use test_suite
  implicit none
  public

  !> Longer than several passes of the longest vector loop the compilers build here (gfortran's integer
  !> loop, 64 values per pass on AVX2), and not a multiple of it, so a remainder loop runs too.
  integer(int32), parameter :: N_VALUES = 269_int32
  integer(int32), parameter :: ARG_POS = 5_int32

contains

  !> Get array of all available tox_errors tests.
  function get_all_tests_tox_errors() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(11))

    all_tests(1) = test_case("test_real_refuses_nan_at_every_position", test_real_refuses_nan_at_every_position)
    all_tests(2) = test_case("test_real_refuses_infinity_at_every_position", test_real_refuses_infinity_at_every_position)
    all_tests(3) = test_case("test_real_out_of_range_at_every_position", test_real_out_of_range_at_every_position)
    all_tests(4) = test_case("test_real_the_first_problem_names_the_error", test_real_the_first_problem_names_the_error)
    all_tests(5) = test_case("test_real_nan_bound_is_ignored_for_finite_values", test_real_nan_bound_is_ignored_for_finite_values)
    all_tests(6) = test_case("test_real_opt_outs_at_every_position", test_real_opt_outs_at_every_position)
    all_tests(7) = test_case("test_real_edge_calls", test_real_edge_calls)
    all_tests(8) = test_case("test_int_out_of_range_at_every_position", test_int_out_of_range_at_every_position)
    all_tests(9) = test_case("test_int_sentinel_at_every_position", test_int_sentinel_at_every_position)
    all_tests(10) = test_case("test_int_edge_calls", test_int_edge_calls)
    all_tests(11) = test_case("test_distance_matrix_refuses_nan_and_infinity", test_distance_matrix_refuses_nan_and_infinity)
  end function get_all_tests_tox_errors

  !> Writes `bad` into a clean array at each position in turn and returns at how many positions the
  !> real validator did not answer `expected`.
  integer(int32) function n_wrong_real(bad, expected, min, max, sentinel, allow_nan, allow_infinite) result(n_wrong)
    real(real64), intent(in) :: bad
    integer(int32), intent(in) :: expected
    real(real64), intent(in), optional :: min, max, sentinel
    logical(c_bool), intent(in), optional :: allow_nan, allow_infinite
    real(real64) :: values(N_VALUES)
    integer(int32) :: i_position, ierr

    n_wrong = 0
    values = 1.0_real64
    do i_position = 1, N_VALUES
      values(i_position) = bad
      ierr = ERR_OK
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, min, max, sentinel, allow_nan, allow_infinite)
      if (ierr /= expected) n_wrong = n_wrong + 1
      values(i_position) = 1.0_real64
    end do
  end function n_wrong_real

  !> The same for the integer validator, on a clean array of ones.
  integer(int32) function n_wrong_int(bad, expected, min, max, sentinel) result(n_wrong)
    integer(int32), intent(in) :: bad, expected
    integer(int32), intent(in), optional :: min, max, sentinel
    integer(int32) :: values(N_VALUES), i_position, ierr

    n_wrong = 0
    values = 1_int32
    do i_position = 1, N_VALUES
      values(i_position) = bad
      ierr = ERR_OK
      call validate_all_in_range_int(values, N_VALUES, ierr, ARG_POS, min, max, sentinel)
      if (ierr /= expected) n_wrong = n_wrong + 1
      values(i_position) = 1_int32
    end do
  end function n_wrong_int

  subroutine test_real_refuses_nan_at_every_position()
    real(real64) :: nan, pos_inf, neg_inf
    integer(int32) :: refused

    nan = ieee_value(1.0_real64, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_real64, ieee_positive_inf)
    neg_inf = ieee_value(1.0_real64, ieee_negative_inf)
    refused = create_err_code(ERR_NAN_INF, ARG_POS)

    call assert_equal_int(n_wrong_real(nan, refused), 0_int32, "NaN, no bounds")
    call assert_equal_int(n_wrong_real(-nan, refused), 0_int32, "-NaN, no bounds")
    call assert_equal_int(n_wrong_real(nan, refused, min=-10.0_real64, max=10.0_real64), 0_int32, "NaN, finite bounds")
    call assert_equal_int(n_wrong_real(nan, refused, min=neg_inf, max=pos_inf), 0_int32, "NaN, infinite bounds")
  end subroutine test_real_refuses_nan_at_every_position

  !> The bounds never let an infinity through, not even infinite ones.
  subroutine test_real_refuses_infinity_at_every_position()
    real(real64) :: pos_inf, neg_inf
    integer(int32) :: refused

    pos_inf = ieee_value(1.0_real64, ieee_positive_inf)
    neg_inf = ieee_value(1.0_real64, ieee_negative_inf)
    refused = create_err_code(ERR_NAN_INF, ARG_POS)

    call assert_equal_int(n_wrong_real(pos_inf, refused), 0_int32, "+Inf, no bounds")
    call assert_equal_int(n_wrong_real(neg_inf, refused), 0_int32, "-Inf, no bounds")
    call assert_equal_int(n_wrong_real(pos_inf, refused, min=-10.0_real64, max=10.0_real64), 0_int32, "+Inf, finite bounds")
    call assert_equal_int(n_wrong_real(neg_inf, refused, min=-10.0_real64, max=10.0_real64), 0_int32, "-Inf, finite bounds")
    call assert_equal_int(n_wrong_real(pos_inf, refused, max=pos_inf), 0_int32, "+Inf, max=+Inf")
    call assert_equal_int(n_wrong_real(neg_inf, refused, min=neg_inf), 0_int32, "-Inf, min=-Inf")
    call assert_equal_int(n_wrong_real(pos_inf, refused, min=neg_inf, max=pos_inf), 0_int32, "+Inf, infinite bounds")
    call assert_equal_int(n_wrong_real(neg_inf, refused, min=neg_inf, max=pos_inf), 0_int32, "-Inf, infinite bounds")
  end subroutine test_real_refuses_infinity_at_every_position

  !> Inclusive bounds: on the bound passes, one float past it does not.
  subroutine test_real_out_of_range_at_every_position()
    real(real64), parameter :: lower = -2.0_real64, upper = 2.0_real64
    real(real64) :: past_lower, past_upper
    integer(int32) :: refused

    past_lower = lower - spacing(lower)
    past_upper = upper + spacing(upper)
    refused = create_err_code(ERR_INVALID_INPUT, ARG_POS)

    call assert_equal_int(n_wrong_real(lower, ERR_OK, min=lower, max=upper), 0_int32, "on the lower bound")
    call assert_equal_int(n_wrong_real(upper, ERR_OK, min=lower, max=upper), 0_int32, "on the upper bound")
    call assert_equal_int(n_wrong_real(past_lower, refused, min=lower, max=upper), 0_int32, "just below the lower bound")
    call assert_equal_int(n_wrong_real(past_upper, refused, min=lower, max=upper), 0_int32, "just above the upper bound")
    call assert_equal_int(n_wrong_real(huge(1.0_real64), ERR_OK), 0_int32, "huge, no bounds")
    call assert_equal_int(n_wrong_real(-huge(1.0_real64), ERR_OK), 0_int32, "-huge, no bounds")
  end subroutine test_real_out_of_range_at_every_position

  !> Of two different problems, the one further to the front decides the error.
  subroutine test_real_the_first_problem_names_the_error()
    real(real64) :: values(N_VALUES), nan
    integer(int32) :: i_first, ierr, n_wrong

    nan = ieee_value(1.0_real64, ieee_quiet_nan)
    n_wrong = 0
    do i_first = 1, N_VALUES - 1
      values = 1.0_real64
      values(i_first) = nan
      values(i_first + 1) = 100.0_real64
      ierr = ERR_OK
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, max=10.0_real64)
      if (ierr /= create_err_code(ERR_NAN_INF, ARG_POS)) n_wrong = n_wrong + 1

      values(i_first) = 100.0_real64
      values(i_first + 1) = nan
      ierr = ERR_OK
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, max=10.0_real64)
      if (ierr /= create_err_code(ERR_INVALID_INPUT, ARG_POS)) n_wrong = n_wrong + 1
    end do
    call assert_equal_int(n_wrong, 0_int32, "the first problem names the error")
  end subroutine test_real_the_first_problem_names_the_error

  !> `val < NaN` is false, so a NaN bound does not limit anything, and NaN values stay refused.
  subroutine test_real_nan_bound_is_ignored_for_finite_values()
    real(real64) :: nan, pos_inf
    integer(int32) :: refused

    nan = ieee_value(1.0_real64, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_real64, ieee_positive_inf)
    refused = create_err_code(ERR_NAN_INF, ARG_POS)

    call assert_equal_int(n_wrong_real(1e300_real64, ERR_OK, min=nan), 0_int32, "finite value, min=NaN")
    call assert_equal_int(n_wrong_real(1e300_real64, ERR_OK, max=nan), 0_int32, "finite value, max=NaN")
    call assert_equal_int(n_wrong_real(1e300_real64, ERR_OK, min=nan, max=nan), 0_int32, "finite value, both NaN")
    call assert_equal_int(n_wrong_real(100.0_real64, create_err_code(ERR_INVALID_INPUT, ARG_POS), min=nan, max=10.0_real64), &
                          0_int32, "the other bound still counts")
    call assert_equal_int(n_wrong_real(nan, refused, min=nan, max=nan), 0_int32, "NaN value, both NaN")
    call assert_equal_int(n_wrong_real(pos_inf, refused, min=nan, max=nan), 0_int32, "+Inf value, both NaN")
  end subroutine test_real_nan_bound_is_ignored_for_finite_values

  !> `sentinel`, `allow_nan` and `allow_infinite` let through what the scan would refuse.
  subroutine test_real_opt_outs_at_every_position()
    real(real64) :: nan, pos_inf, neg_inf
    integer(int32) :: refused

    nan = ieee_value(1.0_real64, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_real64, ieee_positive_inf)
    neg_inf = ieee_value(1.0_real64, ieee_negative_inf)
    refused = create_err_code(ERR_NAN_INF, ARG_POS)

    call assert_equal_int(n_wrong_real(nan, ERR_OK, allow_nan=.true._c_bool), 0_int32, "allow_nan")
    call assert_equal_int(n_wrong_real(nan, refused, allow_nan=.false._c_bool), 0_int32, "allow_nan=.false.")
    call assert_equal_int(n_wrong_real(pos_inf, refused, allow_nan=.true._c_bool), 0_int32, "allow_nan is not allow_infinite")
    call assert_equal_int(n_wrong_real(pos_inf, ERR_OK, allow_infinite=.true._c_bool), 0_int32, "allow_infinite, +Inf")
    call assert_equal_int(n_wrong_real(neg_inf, ERR_OK, allow_infinite=.true._c_bool), 0_int32, "allow_infinite, -Inf")
    call assert_equal_int(n_wrong_real(nan, refused, allow_infinite=.true._c_bool), 0_int32, "allow_infinite is not allow_nan")
    call assert_equal_int(n_wrong_real(nan, ERR_OK, sentinel=nan), 0_int32, "NaN sentinel")
    call assert_equal_int(n_wrong_real(pos_inf, ERR_OK, sentinel=pos_inf), 0_int32, "+Inf sentinel")
    call assert_equal_int(n_wrong_real(-1.0_real64, ERR_OK, min=0.0_real64, sentinel=-1.0_real64), 0_int32, "finite sentinel")
    call assert_equal_int(n_wrong_real(-2.0_real64, create_err_code(ERR_INVALID_INPUT, ARG_POS), min=0.0_real64, &
                                       sentinel=-1.0_real64), 0_int32, "another value than the sentinel")
  end subroutine test_real_opt_outs_at_every_position

  subroutine test_real_edge_calls()
    real(real64) :: values(N_VALUES)
    real(real64), allocatable :: absent_array(:)
    integer(int32) :: ierr, earlier_error

    values = 1.0_real64
    earlier_error = create_err_code(ERR_EMPTY_INPUT, 9_int32)

    ierr = ERR_OK
    call validate_all_in_range_real(values, 0_int32, ierr, ARG_POS, min=5.0_real64, max=-5.0_real64)
    call assert_equal_int(ierr, ERR_OK, "no elements: nothing to refuse")

    ierr = ERR_OK
    call validate_all_in_range_real(absent_array, N_VALUES, ierr, ARG_POS)
    call assert_equal_int(ierr, ERR_OK, "absent array")

    ierr = ERR_OK
    call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, min=5.0_real64, max=-5.0_real64)
    call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, ARG_POS), "min > max refuses everything")

    ierr = earlier_error
    call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS)
    call assert_equal_int(ierr, earlier_error, "clean array: an earlier error stays")

    values(N_VALUES) = ieee_value(1.0_real64, ieee_quiet_nan)
    ierr = earlier_error
    call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS)
    call assert_equal_int(ierr, earlier_error, "bad array: the earlier error stays first")
  end subroutine test_real_edge_calls

  subroutine test_int_out_of_range_at_every_position()
    integer(int32) :: refused

    refused = create_err_code(ERR_INVALID_INPUT, ARG_POS)

    call assert_equal_int(n_wrong_int(-2_int32, ERR_OK, min=-2_int32, max=2_int32), 0_int32, "on the lower bound")
    call assert_equal_int(n_wrong_int(2_int32, ERR_OK, min=-2_int32, max=2_int32), 0_int32, "on the upper bound")
    call assert_equal_int(n_wrong_int(-3_int32, refused, min=-2_int32, max=2_int32), 0_int32, "just below the lower bound")
    call assert_equal_int(n_wrong_int(3_int32, refused, min=-2_int32, max=2_int32), 0_int32, "just above the upper bound")
    call assert_equal_int(n_wrong_int(huge(1_int32), ERR_OK), 0_int32, "huge, no bounds")
    ! one below -huge exists in two's complement, and the default lower bound refuses it
    call assert_equal_int(n_wrong_int(-huge(1_int32) - 1_int32, refused), 0_int32, "most negative integer, no bounds")
  end subroutine test_int_out_of_range_at_every_position

  subroutine test_int_sentinel_at_every_position()
    call assert_equal_int(n_wrong_int(-1_int32, ERR_OK, min=0_int32, sentinel=-1_int32), 0_int32, "the sentinel")
    call assert_equal_int(n_wrong_int(-2_int32, create_err_code(ERR_INVALID_INPUT, ARG_POS), min=0_int32, sentinel=-1_int32), &
                          0_int32, "another value than the sentinel")
  end subroutine test_int_sentinel_at_every_position

  subroutine test_int_edge_calls()
    integer(int32) :: values(N_VALUES), ierr, earlier_error
    integer(int32), allocatable :: absent_array(:)

    values = 1_int32
    earlier_error = create_err_code(ERR_EMPTY_INPUT, 9_int32)

    ierr = ERR_OK
    call validate_all_in_range_int(values, 0_int32, ierr, ARG_POS, min=5_int32, max=-5_int32)
    call assert_equal_int(ierr, ERR_OK, "no elements: nothing to refuse")

    ierr = ERR_OK
    call validate_all_in_range_int(absent_array, N_VALUES, ierr, ARG_POS)
    call assert_equal_int(ierr, ERR_OK, "absent array")

    ierr = ERR_OK
    call validate_all_in_range_int(values, N_VALUES, ierr, ARG_POS, min=5_int32, max=-5_int32)
    call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, ARG_POS), "min > max refuses everything")

    ierr = earlier_error
    call validate_all_in_range_int(values, N_VALUES, ierr, ARG_POS)
    call assert_equal_int(ierr, earlier_error, "clean array: an earlier error stays")

    values(N_VALUES) = 10_int32
    ierr = earlier_error
    call validate_all_in_range_int(values, N_VALUES, ierr, ARG_POS, max=5_int32)
    call assert_equal_int(ierr, earlier_error, "bad array: the earlier error stays first")
  end subroutine test_int_edge_calls

  !> The distance matrix validator reaches the scan through the real one.
  subroutine test_distance_matrix_refuses_nan_and_infinity()
    integer(int32), parameter :: n = 12_int32
    real(real64) :: distances(n, n), bad_values(2)
    integer(int32) :: i_bad, i_position, ierr, n_wrong

    bad_values = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]

    distances = 1.0_real64
    do i_position = 1, n
      distances(i_position, i_position) = 0.0_real64
    end do
    ierr = ERR_OK
    call validate_distance_matrix(distances, n, ierr, ARG_POS)
    call assert_equal_int(ierr, ERR_OK, "clean matrix")

    n_wrong = 0
    do i_bad = 1, size(bad_values)
      do i_position = 1, n*n
        block
          real(real64) :: flat(n*n)
          flat = reshape(distances, [n*n])
          flat(i_position) = bad_values(i_bad)
          ierr = ERR_OK
          call validate_distance_matrix(reshape(flat, [n, n]), n, ierr, ARG_POS)
          if (ierr /= create_err_code(ERR_NAN_INF, ARG_POS)) n_wrong = n_wrong + 1
        end block
      end do
    end do
    call assert_equal_int(n_wrong, 0_int32, "NaN and +Inf at every position of the matrix")
  end subroutine test_distance_matrix_refuses_nan_and_infinity

end module mod_test_tox_errors
