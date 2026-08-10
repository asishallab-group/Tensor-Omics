!> Unit test suite for tox_shape_truthful_clustering_accept (accept_ensemble), generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_accept_kernel.F90.
module mod_test_shape_truthful_clustering_accept
  use tox_shape_truthful_clustering_accept, only: accept_ensemble_alloc
  use tox_errors, only: is_ok, is_err
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use test_suite, only: test_case
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_shape_truthful_clustering_accept() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate (all_tests(17))

    all_tests(1) = test_case("test_accept_ensemble_identical", test_accept_ensemble_identical)
    all_tests(2) = test_case("test_accept_ensemble_chordal_exceeds_max", test_accept_ensemble_chordal_exceeds_max)
    all_tests(3) = test_case("test_accept_ensemble_chordal_within_max", test_accept_ensemble_chordal_within_max)
    all_tests(4) = test_case("test_accept_ensemble_rejects_cumulative_drift_from_first", &
                              test_accept_ensemble_rejects_cumulative_drift_from_first)
    all_tests(5) = test_case("test_accept_ensemble_accepts_small_drift_from_both", &
                              test_accept_ensemble_accepts_small_drift_from_both)
    all_tests(6) = test_case("test_accept_ensemble_d_to_first_exceeds_dmax", &
                              test_accept_ensemble_d_to_first_exceeds_dmax)
    all_tests(7) = test_case("test_accept_ensemble_d_two_fold_within_dmax", &
                              test_accept_ensemble_d_two_fold_within_dmax)
    all_tests(8) = test_case("test_accept_ensemble_g_ratio_exceeds_max", test_accept_ensemble_g_ratio_exceeds_max)
    all_tests(9) = test_case("test_accept_ensemble_g_ratio_within_max", test_accept_ensemble_g_ratio_within_max)
    all_tests(10) = test_case("test_accept_ensemble_rmse_ratio_exceeds_max", &
                               test_accept_ensemble_rmse_ratio_exceeds_max)
    all_tests(11) = test_case("test_accept_ensemble_rmse_ratio_within_max", &
                               test_accept_ensemble_rmse_ratio_within_max)
    all_tests(12) = test_case("test_accept_ensemble_nonpositive_g", test_accept_ensemble_nonpositive_g)
    all_tests(13) = test_case("test_accept_ensemble_nonpositive_normal_error", &
                               test_accept_ensemble_nonpositive_normal_error)
    all_tests(14) = test_case("test_accept_ensemble_zero_normal_error_is_accepted", &
                               test_accept_ensemble_zero_normal_error_is_accepted)
    all_tests(15) = test_case("test_accept_ensemble_d_first_out_of_range", test_accept_ensemble_d_first_out_of_range)
    all_tests(16) = test_case("test_accept_ensemble_chordal_frac_out_of_range", &
                               test_accept_ensemble_chordal_frac_out_of_range)
    all_tests(17) = test_case("test_accept_ensemble_history_len_out_of_range", &
                               test_accept_ensemble_history_len_out_of_range)
  end function get_all_tests_shape_truthful_clustering_accept

  !> Identical bases at every retained iteration (U_first, U_history, U_tp1), identical d, G,
  !| normal_error: all four criteria trivially satisfied, even at zero tolerance for d, G, RMSE.
  subroutine test_accept_ensemble_identical()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               5.0d0, 1.0d0, U, 2_int32, 5.0d0, 1.0d0, &
                               0.1d0, 0_int32, 0.0d0, 0.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: identical observables must be accepted")
  end subroutine test_accept_ensemble_identical

  !> A 60-degree rotation of a 1D tangent basis (U_first = U_history(1), so this isolates the
  !| chordal-distance formula/threshold itself from the cumulative-drift machinery): chordal
  !| distance = sin(60deg) ~ 0.866, rejected at a fraction-of-range of 0.5.
  subroutine test_accept_ensemble_chordal_exceeds_max()
    integer(int32), parameter :: d_dim = 2, o = 1
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    real(real64)   :: pi, theta
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0   ! 60 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    U_history(:, :, 1) = U_t
    d_history(1)       = 1_int32

    call accept_ensemble_alloc(d_dim, o, U_t, 1_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U_tp1, 1_int32, 1.0d0, 1.0d0, &
                               0.5d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, &
      "accept_ensemble: a 60-degree tangent rotation must be rejected at chordal_dist_max_as_prcnt_of_range=0.5")
  end subroutine test_accept_ensemble_chordal_exceeds_max

  !> The same 60-degree rotation, accepted once the fraction-of-range is raised to 0.9
  !| (sin(60deg) ~ 0.866 <= 0.9).
  subroutine test_accept_ensemble_chordal_within_max()
    integer(int32), parameter :: d_dim = 2, o = 1
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    real(real64)   :: pi, theta
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0        ! 60 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    U_history(:, :, 1) = U_t
    d_history(1)       = 1_int32

    call accept_ensemble_alloc(d_dim, o, U_t, 1_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U_tp1, 1_int32, 1.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, &
      "accept_ensemble: a 60-degree tangent rotation must be accepted at chordal_dist_max_as_prcnt_of_range=0.9")
  end subroutine test_accept_ensemble_chordal_within_max

  !> The P5 regression this whole redesign fixes: a candidate only 5 degrees from the most
  !| recently accepted state (U_history(1), at 75deg) -- which a step-to-step-only check would
  !| accept -- but 80 degrees from the ensemble's own bootstrap state (U_first, at 0deg).
  !| Comparing against the reference set (not just the previous state) must reject this.
  subroutine test_accept_ensemble_rejects_cumulative_drift_from_first()
    integer(int32), parameter :: d_dim = 2, o = 1
    real(real64)   :: U_first(d_dim, d_dim), U_tp1(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    real(real64)   :: pi
    logical        :: is_accepted
    integer(int32) :: ierr

    pi = 4.0d0 * atan(1.0d0)

    U_first = 0.0d0
    U_first(1,1) = 1.0d0
    U_first(2,2) = 1.0d0

    U_history(:, :, 1)   = 0.0d0
    U_history(1,1,1)     = cos(75.0d0 * pi / 180.0d0)
    U_history(2,1,1)     = sin(75.0d0 * pi / 180.0d0)
    U_history(2,2,1)     = 1.0d0
    d_history(1)         = 1_int32

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(80.0d0 * pi / 180.0d0)
    U_tp1(2,1) = sin(80.0d0 * pi / 180.0d0)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, o, U_first, 1_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U_tp1, 1_int32, 1.0d0, 1.0d0, &
                               0.5d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, &
      "accept_ensemble: a candidate close to the previous state but far from U_first must be rejected")
  end subroutine test_accept_ensemble_rejects_cumulative_drift_from_first

  !> Companion to the above: small drift from both U_first (6deg) and U_history(1) (3deg) is
  !| accepted at the same threshold.
  subroutine test_accept_ensemble_accepts_small_drift_from_both()
    integer(int32), parameter :: d_dim = 2, o = 1
    real(real64)   :: U_first(d_dim, d_dim), U_tp1(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    real(real64)   :: pi
    logical        :: is_accepted
    integer(int32) :: ierr

    pi = 4.0d0 * atan(1.0d0)

    U_first = 0.0d0
    U_first(1,1) = 1.0d0
    U_first(2,2) = 1.0d0

    U_history(:, :, 1)   = 0.0d0
    U_history(1,1,1)     = cos(3.0d0 * pi / 180.0d0)
    U_history(2,1,1)     = sin(3.0d0 * pi / 180.0d0)
    U_history(2,2,1)     = 1.0d0
    d_history(1)         = 1_int32

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(6.0d0 * pi / 180.0d0)
    U_tp1(2,1) = sin(6.0d0 * pi / 180.0d0)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, o, U_first, 1_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U_tp1, 1_int32, 1.0d0, 1.0d0, &
                               0.5d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, &
      "accept_ensemble: small drift from both U_first and U_history(1) must be accepted")
  end subroutine test_accept_ensemble_accepts_small_drift_from_both

  !> d_to_last = |d_tp1 - d_history(1)| = 0 (fine on its own), but d_to_first =
  !| |d_tp1 - d_first| = 2 > d_max=1: a d_to_last-only check would wrongly accept this: the
  !| two-fold max() must reject it. d_first=0 makes the chordal criterion vacuous against it
  !| (no shared dimension), so only criterion (2) is at play here.
  subroutine test_accept_ensemble_d_to_first_exceeds_dmax()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 0_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 1_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, &
      "accept_ensemble: max(d_to_first, d_to_last) must reject even when d_to_last alone is fine")
  end subroutine test_accept_ensemble_d_to_first_exceeds_dmax

  !> Same setup, accepted once d_max is raised to 2 (>= max(d_to_first, d_to_last) = 2).
  subroutine test_accept_ensemble_d_two_fold_within_dmax()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 0_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 2_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: d_max=2 must tolerate max(d_to_first, d_to_last)=2")
  end subroutine test_accept_ensemble_d_two_fold_within_dmax

  !> G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0. Consecutive-only
  !| (against G_t, the most recently accepted iteration), unaffected by history.
  subroutine test_accept_ensemble_g_ratio_exceeds_max()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 10.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d0, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a 10x change in G must be rejected at G_max=1.0")
  end subroutine test_accept_ensemble_g_ratio_exceeds_max

  !> Same 10x change in G, accepted once G_max is raised past ln(10) ~ 2.303.
  subroutine test_accept_ensemble_g_ratio_within_max()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 10.0d0, 1.0d0, &
                               0.9d0, 0_int32, 3.0d0, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: a 10x change in G must be accepted at G_max=3.0")
  end subroutine test_accept_ensemble_g_ratio_within_max

  !> normal_error changes by a factor of 10 (RMSE by sqrt(10) ~ 3.162, log ~ 1.151): rejected
  !| at RMSE_change_max=1.0.
  subroutine test_accept_ensemble_rmse_ratio_exceeds_max()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 10.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, &
      "accept_ensemble: a 10x change in normal_error must be rejected at RMSE_change_max=1.0")
  end subroutine test_accept_ensemble_rmse_ratio_exceeds_max

  !> Same 10x change in normal_error, accepted once RMSE_change_max is raised past
  !| |log(sqrt(10))| ~ 1.151.
  subroutine test_accept_ensemble_rmse_ratio_within_max()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 10.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.2d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, &
      "accept_ensemble: a 10x change in normal_error must be accepted at RMSE_change_max=1.2")
  end subroutine test_accept_ensemble_rmse_ratio_within_max

  !> Non-positive spectral gaps make log(G_tp1/G_t) undefined and must be rejected by validation.
  subroutine test_accept_ensemble_nonpositive_g()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               0.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject G_t <= 0")
  end subroutine test_accept_ensemble_nonpositive_g

  !> Negative normal_error is physically impossible (a sum of eigenvalues) and must be
  !| rejected by validation -- unlike G_t, zero itself is valid (see
  !| test_accept_ensemble_zero_normal_error_is_accepted below), since normal_error, unlike G,
  !| carries no ratio-internal +epsilon protection of its own.
  subroutine test_accept_ensemble_nonpositive_normal_error()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, -1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject normal_error_t < 0")
  end subroutine test_accept_ensemble_nonpositive_normal_error

  !> Zero normal_error at both t and t+1 (a perfectly flat/collinear ensemble, no noise at all
  !| in the normal directions) must not crash or spuriously reject via a log(0/0) -- the
  !| +epsilon(1.0_real64) guard inside the RMSE ratio (mirroring observable's own G(r)
  !| epsilon convention) keeps it well-defined and accepted.
  subroutine test_accept_ensemble_zero_normal_error_is_accepted()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 0.0d0, U, 2_int32, 1.0d0, 0.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: zero normal_error at both t and t+1 must be accepted, not NaN-rejected")
  end subroutine test_accept_ensemble_zero_normal_error_is_accepted

  !> d_first out of [0, n_dimensions] must be rejected by validation.
  subroutine test_accept_ensemble_d_first_out_of_range()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, d_dim + 1_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject d_first > n_dimensions")
  end subroutine test_accept_ensemble_d_first_out_of_range

  !> chordal_dist_max_as_prcnt_of_range out of [0, 1] must be rejected by validation.
  subroutine test_accept_ensemble_chordal_frac_out_of_range()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 1_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               1.5d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject chordal_dist_max_as_prcnt_of_range > 1")
  end subroutine test_accept_ensemble_chordal_frac_out_of_range

  !> history_len out of [1, o] must be rejected by validation.
  subroutine test_accept_ensemble_history_len_out_of_range()
    integer(int32), parameter :: d_dim = 3, o = 1
    real(real64)   :: U(d_dim, d_dim), U_history(d_dim, d_dim, o)
    integer(int32) :: d_history(o)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U = 0.0d0
    do i = 1, d_dim
      U(i, i) = 1.0d0
    end do
    U_history(:, :, 1) = U
    d_history(1)       = 2_int32

    call accept_ensemble_alloc(d_dim, o, U, 2_int32, U_history, d_history, 2_int32, &
                               1.0d0, 1.0d0, U, 2_int32, 1.0d0, 1.0d0, &
                               0.9d0, 0_int32, 1.0d10, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject history_len > o")
  end subroutine test_accept_ensemble_history_len_out_of_range

end module mod_test_shape_truthful_clustering_accept
