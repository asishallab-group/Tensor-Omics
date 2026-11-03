!> Unit test suite for the array validators in tox_errors.
!> They answer from a fast scan when an array is clean and from a per-element loop otherwise, so every
!> test here puts the problem at every position of an array long enough to span several vector
!> lengths: the vector body, its remainder and the seam between them must all see it. This is also
!> the guard against a compiler assuming NaN cannot occur, so it matters under every build profile.
!> Every check also requires that no IEEE flag was raised: a NaN must stay quiet, so a host that
!> halts on invalid gets `ierr` back instead of stopping inside a validator (#218).

module mod_test_tox_errors
  use asserts
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use, intrinsic :: ieee_exceptions, only: ieee_all, ieee_get_flag, ieee_set_flag
  use test_suite
  implicit none
  public

  !> Longer than several passes of the longest vector loop the compilers build here (gfortran's integer
  !> loop, 64 values per pass on AVX2), and not a multiple of it, so a remainder loop runs too.
  integer(int32), parameter :: N_VALUES = 269_int32
  integer(int32), parameter :: ARG_POS = 5_int32
  !> Built from its bits: `ieee_value` may hand back a quietened one.
  real(real64), parameter :: SIGNALLING_NAN = transfer(int(z'7FF0000000000001', int64), 1.0_real64)
  !> The same with the sign bit set (gfortran refuses to negate a NaN constant).
  real(real64), parameter :: NEGATIVE_SIGNALLING_NAN = &
    transfer(ior(int(z'7FF0000000000001', int64), ibset(0_int64, 63)), 1.0_real64)

contains

  !> Get array of all available tox_errors tests.
  function get_all_tests_tox_errors() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(14))

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
    all_tests(12) = test_case("test_real_refuses_signalling_nan_at_every_position", &
                              test_real_refuses_signalling_nan_at_every_position)
    all_tests(13) = test_case("test_real_nan_bounds_raise_no_flag", test_real_nan_bounds_raise_no_flag)
    all_tests(14) = test_case("test_real_scalar_raises_no_flag", test_real_scalar_raises_no_flag)
  end function get_all_tests_tox_errors

  !> Whether an IEEE flag was raised since the flags were last cleared; clears them for the next check.
  logical function raised_a_flag()
    logical :: raised(size(ieee_all))

    call ieee_get_flag(ieee_all, raised)
    raised_a_flag = any(raised)
    call ieee_set_flag(ieee_all, .false.)
  end function raised_a_flag

  !> Writes `bad` into a clean array at each position in turn and returns at how many positions the
  !> real validator did not answer `expected`, or raised an IEEE flag.
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
      call ieee_set_flag(ieee_all, .false.)
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, min, max, sentinel, allow_nan, allow_infinite)
      if (ierr /= expected .or. raised_a_flag()) n_wrong = n_wrong + 1
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
    ! -0.0 == 0.0, so each is on the other as a lower bound (the array's other values are ones)
    call assert_equal_int(n_wrong_real(sign(0.0_real64, -1.0_real64), ERR_OK, min=0.0_real64), 0_int32, "-0.0 on min=0.0")
    call assert_equal_int(n_wrong_real(0.0_real64, ERR_OK, min=sign(0.0_real64, -1.0_real64)), 0_int32, "0.0 on min=-0.0")
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
      call ieee_set_flag(ieee_all, .false.)
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, max=10.0_real64)
      if (ierr /= create_err_code(ERR_NAN_INF, ARG_POS) .or. raised_a_flag()) n_wrong = n_wrong + 1

      values(i_first) = 100.0_real64
      values(i_first + 1) = nan
      ierr = ERR_OK
      call ieee_set_flag(ieee_all, .false.)
      call validate_all_in_range_real(values, N_VALUES, ierr, ARG_POS, max=10.0_real64)
      if (ierr /= create_err_code(ERR_INVALID_INPUT, ARG_POS) .or. raised_a_flag()) n_wrong = n_wrong + 1
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
    real(real64) :: distances(n, n), bad_values(3)
    integer(int32) :: i_bad, i_position, ierr, n_wrong

    bad_values = [ieee_value(1.0_real64, ieee_quiet_nan), SIGNALLING_NAN, ieee_value(1.0_real64, ieee_positive_inf)]

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
          call ieee_set_flag(ieee_all, .false.)
          call validate_distance_matrix(reshape(flat, [n, n]), n, ierr, ARG_POS)
          if (ierr /= create_err_code(ERR_NAN_INF, ARG_POS) .or. raised_a_flag()) n_wrong = n_wrong + 1
        end block
      end do
    end do
    call assert_equal_int(n_wrong, 0_int32, "NaN, signalling NaN and +Inf at every position of the matrix")
  end subroutine test_distance_matrix_refuses_nan_and_infinity

  !> A signalling NaN is refused like a quiet one, and raises no flag either: every comparison of reals
  !> would raise invalid on it, and on gfortran so would ieee_is_nan.
  subroutine test_real_refuses_signalling_nan_at_every_position()
    integer(int32) :: refused

    refused = create_err_code(ERR_NAN_INF, ARG_POS)

    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, refused), 0_int32, "sNaN, no bounds")
    call assert_equal_int(n_wrong_real(NEGATIVE_SIGNALLING_NAN, refused), 0_int32, "-sNaN, no bounds")
    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, refused, min=-10.0_real64, max=10.0_real64), 0_int32, &
                          "sNaN, finite bounds")
    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, ERR_OK, allow_nan=.true._c_bool), 0_int32, "sNaN, allow_nan")
    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, refused, allow_infinite=.true._c_bool), 0_int32, &
                          "sNaN, allow_infinite")
    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, ERR_OK, sentinel=SIGNALLING_NAN), 0_int32, "sNaN sentinel")
    call assert_equal_int(n_wrong_real(SIGNALLING_NAN, ERR_OK, sentinel=ieee_value(1.0_real64, ieee_quiet_nan)), 0_int32, &
                          "sNaN value, quiet NaN sentinel")
    call assert_equal_int(n_wrong_real(2.0_real64, create_err_code(ERR_INVALID_INPUT, ARG_POS), &
                                       max=1.5_real64, sentinel=SIGNALLING_NAN), 0_int32, "finite value, sNaN sentinel")
  end subroutine test_real_refuses_signalling_nan_at_every_position

  !> A NaN bound, quiet or signalling, limits nothing and raises no flag, on the fast path and the per-element one.
  subroutine test_real_nan_bounds_raise_no_flag()
    real(real64) :: nan, pos_inf

    nan = ieee_value(1.0_real64, ieee_quiet_nan)
    pos_inf = ieee_value(1.0_real64, ieee_positive_inf)

    call assert_equal_int(n_wrong_real(1e300_real64, ERR_OK, min=SIGNALLING_NAN), 0_int32, "finite value, min=sNaN")
    call assert_equal_int(n_wrong_real(-1e300_real64, ERR_OK, max=SIGNALLING_NAN), 0_int32, "finite value, max=sNaN")
    call assert_equal_int(n_wrong_real(-1e300_real64, ERR_OK, min=nan, max=SIGNALLING_NAN), 0_int32, "finite value, NaN bounds")
    call assert_equal_int(n_wrong_real(pos_inf, create_err_code(ERR_NAN_INF, ARG_POS), min=SIGNALLING_NAN, max=nan), 0_int32, &
                          "+Inf value, NaN bounds")
    call assert_equal_int(n_wrong_real(-1e300_real64, ERR_OK, min=nan, allow_nan=.true._c_bool), 0_int32, &
                          "finite value, NaN bound, per-element path")
  end subroutine test_real_nan_bounds_raise_no_flag

  !> The scalar validator, which the array one falls back on: no flag for any NaN or infinity, as value, bound or sentinel.
  subroutine test_real_scalar_raises_no_flag()
    real(real64) :: specials(6)
    integer(int32) :: i_value, i_other, ierr, n_loud

    specials = [ieee_value(1.0_real64, ieee_quiet_nan), -ieee_value(1.0_real64, ieee_quiet_nan), SIGNALLING_NAN, &
                NEGATIVE_SIGNALLING_NAN, ieee_value(1.0_real64, ieee_positive_inf), ieee_value(1.0_real64, ieee_negative_inf)]
    n_loud = 0
    call ieee_set_flag(ieee_all, .false.)
    do i_value = 1, size(specials)
      do i_other = 1, size(specials)
        ierr = ERR_OK
        call validate_in_range_real(specials(i_value), ierr, ARG_POS, min=specials(i_other))
        call validate_in_range_real(specials(i_value), ierr, ARG_POS, max=specials(i_other))
        call validate_in_range_real(specials(i_value), ierr, ARG_POS, sentinel=specials(i_other))
        call validate_in_range_real(1.0_real64, ierr, ARG_POS, min=specials(i_other), max=specials(i_value))
        call validate_in_range_real(1.0_real64, ierr, ARG_POS, sentinel=specials(i_other))
        call validate_in_range_real(specials(i_value), ierr, ARG_POS, allow_nan=.true._c_bool, &
                                    allow_infinite=.true._c_bool, sentinel=specials(i_other))
        if (raised_a_flag()) n_loud = n_loud + 1
      end do
    end do
    call assert_equal_int(n_loud, 0_int32, "combinations that raised an IEEE flag")
  end subroutine test_real_scalar_raises_no_flag

end module mod_test_tox_errors
