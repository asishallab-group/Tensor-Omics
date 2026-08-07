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
    allocate (all_tests(8))

    all_tests(1) = test_case("test_accept_ensemble_identical", test_accept_ensemble_identical)
    all_tests(2) = test_case("test_accept_ensemble_angle_exceeds_max", test_accept_ensemble_angle_exceeds_max)
    all_tests(3) = test_case("test_accept_ensemble_angle_within_max", test_accept_ensemble_angle_within_max)
    all_tests(4) = test_case("test_accept_ensemble_d_mismatch_within_dmax", test_accept_ensemble_d_mismatch_within_dmax)
    all_tests(5) = test_case("test_accept_ensemble_d_mismatch_exceeds_dmax", test_accept_ensemble_d_mismatch_exceeds_dmax)
    all_tests(6) = test_case("test_accept_ensemble_g_ratio_exceeds_max", test_accept_ensemble_g_ratio_exceeds_max)
    all_tests(7) = test_case("test_accept_ensemble_nonpositive_g", test_accept_ensemble_nonpositive_g)
    all_tests(8) = test_case("test_accept_ensemble_d_out_of_range", test_accept_ensemble_d_out_of_range)
  end function get_all_tests_shape_truthful_clustering_accept

  !> Identical basis, identical d, identical G: all three criteria trivially
  !| satisfied, even at zero tolerance for d and G.
  subroutine test_accept_ensemble_identical()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 5.0d0, U_tp1, 2_int32, 5.0d0, &
                         0.1d0, 0_int32, 0.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: identical observables must be accepted")
  end subroutine test_accept_ensemble_identical

  !> A 60-degree rotation of a 1D tangent basis, rejected by an alpha_max of 30 degrees.
  subroutine test_accept_ensemble_angle_exceeds_max()
    integer(int32), parameter :: d_dim = 2
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), pi, theta, alpha_max
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0   ! 60 degrees
    alpha_max = pi / 6.0d0  ! 30 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, U_t, 1_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                         alpha_max, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a 60-degree tangent rotation must be rejected at alpha_max=30deg")
  end subroutine test_accept_ensemble_angle_exceeds_max

  !> The same 60-degree rotation, accepted once alpha_max is raised to 70 degrees.
  subroutine test_accept_ensemble_angle_within_max()
    integer(int32), parameter :: d_dim = 2
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), pi, theta, alpha_max
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0        ! 60 degrees
    alpha_max = 7.0d0 * pi / 18.0d0  ! 70 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, U_t, 1_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                         alpha_max, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: a 60-degree tangent rotation must be accepted at alpha_max=70deg")
  end subroutine test_accept_ensemble_angle_within_max

  !> d changed by 1 (2 -> 1): the angle criterion is vacuously satisfied
  !| (no common dimension to compare), and d_max=1 tolerates the change.
  subroutine test_accept_ensemble_d_mismatch_within_dmax()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                         0.1d0, 1_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: a change in d of 1 must be accepted when d_max=1")
  end subroutine test_accept_ensemble_d_mismatch_within_dmax

  !> The same change in d (2 -> 1), rejected once d_max=0.
  subroutine test_accept_ensemble_d_mismatch_exceeds_dmax()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                         0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a change in d of 1 must be rejected when d_max=0")
  end subroutine test_accept_ensemble_d_mismatch_exceeds_dmax

  !> G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
  subroutine test_accept_ensemble_g_ratio_exceeds_max()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 2_int32, 10.0d0, &
                         0.1d0, 0_int32, 1.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a 10x change in G must be rejected at G_max=1.0")
  end subroutine test_accept_ensemble_g_ratio_exceeds_max

  !> Non-positive spectral gaps make log(G_tp1/G_t) undefined and must be rejected by validation.
  subroutine test_accept_ensemble_nonpositive_g()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 0.0d0, U_tp1, 2_int32, 1.0d0, &
                         0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject G_t <= 0")
  end subroutine test_accept_ensemble_nonpositive_g

  !> d_t out of [0, n_dimensions] must be rejected by validation.
  subroutine test_accept_ensemble_d_out_of_range()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, d_dim + 1, 1.0d0, U_tp1, 2_int32, 1.0d0, &
                         0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject d_t > n_dimensions")
  end subroutine test_accept_ensemble_d_out_of_range

end module mod_test_shape_truthful_clustering_accept
