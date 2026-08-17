!> Unit test suite for tox_shape_truthful_clustering_observable (normal_error, tangent_scales,
!| observable), generated from
!| src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_observable_impl.F90.
module mod_test_shape_truthful_clustering_observable
  use tox_shape_truthful_clustering_observable, only: normal_error, tangent_scales, observable, &
    ensemble_final_observable
  use tox_errors, only: is_ok, is_err
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
  use test_suite, only: test_case
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_shape_truthful_clustering_observable() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate (all_tests(20))

    all_tests(1) = test_case("test_normal_error_basic", test_normal_error_basic)
    all_tests(2) = test_case("test_normal_error_zero_tangent_dims", test_normal_error_zero_tangent_dims)
    all_tests(3) = test_case("test_normal_error_all_tangent_dims", test_normal_error_all_tangent_dims)
    all_tests(4) = test_case("test_normal_error_d_out_of_range", test_normal_error_d_out_of_range)
    all_tests(5) = test_case("test_normal_error_negative_eigenvalue", test_normal_error_negative_eigenvalue)
    all_tests(6) = test_case("test_normal_error_zero_dimensions", test_normal_error_zero_dimensions)
    all_tests(7) = test_case("test_tangent_scales_basic", test_tangent_scales_basic)
    all_tests(8) = test_case("test_tangent_scales_all_dims", test_tangent_scales_all_dims)
    all_tests(9) = test_case("test_tangent_scales_zero_dims", test_tangent_scales_zero_dims)
    all_tests(10) = test_case("test_tangent_scales_d_out_of_range", test_tangent_scales_d_out_of_range)
    all_tests(11) = test_case("test_tangent_scales_negative_eigenvalue", test_tangent_scales_negative_eigenvalue)
    all_tests(12) = test_case("test_observable_full_rank_rectangle", test_observable_full_rank_rectangle)
    all_tests(13) = test_case("test_observable_low_rank_padding", test_observable_low_rank_padding)
    all_tests(14) = test_case("test_observable_too_few_members", test_observable_too_few_members)
    all_tests(15) = test_case("test_observable_dimension_too_small", test_observable_dimension_too_small)
    all_tests(16) = test_case("test_ensemble_final_observable_trailing_rejected_column", &
                              test_ensemble_final_observable_trailing_rejected_column)
    all_tests(17) = test_case("test_ensemble_final_observable_no_rejection", &
                              test_ensemble_final_observable_no_rejection)
    all_tests(18) = test_case("test_ensemble_final_observable_has_final_false", &
                              test_ensemble_final_observable_has_final_false)
    all_tests(19) = test_case("test_ensemble_final_observable_multi_ensemble_independence", &
                              test_ensemble_final_observable_multi_ensemble_independence)
    all_tests(20) = test_case("test_ensemble_final_observable_small_o_evicts_accepted", &
                              test_ensemble_final_observable_small_o_evicts_accepted)
  end function get_all_tests_shape_truthful_clustering_observable

  ! --- normal_error -----------------------------------------------------

  !> Basic case: D=3, d=1 -- normal_error is the sum of the two smallest
  !| (normal-space) eigenvalues.
  subroutine test_normal_error_basic()
    integer(int32), parameter :: d_dim = 3, d = 1
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 5.0d0, 1.0d-12, "normal_error basic case (d=1)")
  end subroutine test_normal_error_basic

  !> d=0: no tangent directions, everything is "normal" -- the full sum.
  subroutine test_normal_error_zero_tangent_dims()
    integer(int32), parameter :: d_dim = 3, d = 0
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 14.0d0, 1.0d-12, "normal_error with d=0 (full sum)")
  end subroutine test_normal_error_zero_tangent_dims

  !> d=D: every direction is tangent, nothing left over -- sum over the
  !| empty range must be exactly zero.
  subroutine test_normal_error_all_tangent_dims()
    integer(int32), parameter :: d_dim = 3, d = 3
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 0.0d0, 1.0d-12, "normal_error with d=D (empty sum)")
  end subroutine test_normal_error_all_tangent_dims

  !> d > D must be rejected by validation, not silently read out of bounds.
  subroutine test_normal_error_d_out_of_range()
    integer(int32), parameter :: d_dim = 3, d = 4
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject d > n_dimensions")
  end subroutine test_normal_error_d_out_of_range

  !> A negative eigenvalue is not a valid covariance eigenvalue and must be
  !| rejected by validation.
  subroutine test_normal_error_negative_eigenvalue()
    integer(int32), parameter :: d_dim = 3, d = 1
    real(real64) :: eigenvalues(d_dim) = [9.0d0, -1.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject a negative eigenvalue")
  end subroutine test_normal_error_negative_eigenvalue

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_normal_error_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, d = 0
    real(real64) :: eigenvalues(d_dim)
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject n_dimensions=0")
  end subroutine test_normal_error_zero_dimensions

  ! --- tangent_scales ----------------------------------------------------

  !> Basic case: D=3, d=2 -- tangent_scales is the square root of the two
  !| largest (tangent-space) eigenvalues.
  subroutine test_tangent_scales_basic()
    integer(int32), parameter :: d_dim = 3, d = 2
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    real(real64) :: expected(d) = [3.0d0, 2.0d0]
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(tangent_scales_value, expected, d, 1.0d-12, "tangent_scales basic case (d=2)")
  end subroutine test_tangent_scales_basic

  !> d=D: every direction is tangent.
  subroutine test_tangent_scales_all_dims()
    integer(int32), parameter :: d_dim = 3, d = 3
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    real(real64) :: expected(d) = [3.0d0, 2.0d0, 1.0d0]
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(tangent_scales_value, expected, d, 1.0d-12, "tangent_scales with d=D")
  end subroutine test_tangent_scales_all_dims

  !> d=0: no tangent directions -- must return a well-defined, empty (size
  !| zero) array rather than fail.
  subroutine test_tangent_scales_zero_dims()
    integer(int32), parameter :: d_dim = 3, d = 0
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_int(size(tangent_scales_value), 0, "tangent_scales with d=0 should be an empty array")
  end subroutine test_tangent_scales_zero_dims

  !> d > D must be rejected by validation, not silently read out of bounds.
  subroutine test_tangent_scales_d_out_of_range()
    integer(int32), parameter :: d_dim = 3, d = 4
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "tangent_scales should reject d > n_dimensions")
  end subroutine test_tangent_scales_d_out_of_range

  !> A negative eigenvalue is not a valid covariance eigenvalue and must be
  !| rejected by validation.
  subroutine test_tangent_scales_negative_eigenvalue()
    integer(int32), parameter :: d_dim = 3, d = 2
    real(real64) :: eigenvalues(d_dim) = [9.0d0, -4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "tangent_scales should reject a negative eigenvalue")
  end subroutine test_tangent_scales_negative_eigenvalue

  ! --- observable ------------------------------------------------------

  !> A rectangle in the z=0 plane, embedded in 3D: full economy-mode rank
  !| (rank = min(D,k) = 3 = D, no zero-padding), with distinct, hand-
  !| computable eigenvalues -- the x/y cross term is exactly zero for a
  !| centered rectangle, so U is deterministically the standard basis
  !| (up to a per-column sign flip, which SVD never fixes).
  subroutine test_observable_full_rank_rectangle()
    integer(int32), parameter :: d_dim = 3, n = 4
    real(real64)   :: vectors(d_dim, n)
    logical(c_bool)        :: member_selection_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr
    real(real64), parameter :: tol = 1.0d-6

    vectors(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    vectors(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    vectors(:, 3) = [0.0d0, 1.0d0, 0.0d0]
    vectors(:, 4) = [2.0d0, 1.0d0, 0.0d0]
    member_selection_mask = .true.

    call observable(vectors, d_dim, n, member_selection_mask, n, &
                    U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(mu, [1.0d0, 0.5d0, 0.0d0], d_dim, tol, "observable: rectangle mean")
    call assert_equal_int(d, 2, "observable: rectangle intrinsic dimension")
    call assert_equal_real(eigenvalues(1), 4.0d0/3.0d0, tol, "observable: rectangle lambda_1")
    call assert_equal_real(eigenvalues(2), 1.0d0/3.0d0, tol, "observable: rectangle lambda_2")
    call assert_equal_real(eigenvalues(3), 0.0d0, tol, "observable: rectangle lambda_3")
    call assert_equal_real(normal_error_value, 0.0d0, tol, "observable: rectangle normal_error")
    call assert_equal_real(tangent_scales_value(1), sqrt(4.0d0/3.0d0), tol, "observable: rectangle tangent_scales(1)")
    call assert_equal_real(tangent_scales_value(2), sqrt(1.0d0/3.0d0), tol, "observable: rectangle tangent_scales(2)")
    call assert_equal_real(tangent_scales_value(3), 0.0d0, tol, "observable: rectangle tangent_scales(3)")
    call assert_true(G > 1.0d10, "observable: rectangle spectral gap should be huge at the true rank boundary")

    ! U columns are the standard basis up to sign (diagonal covariance with
    ! distinct eigenvalues 4, 1, 0 leaves no rotational freedom).
    call assert_equal_real(abs(U(1,1)), 1.0d0, tol, "observable: rectangle U column 1, x")
    call assert_equal_real(abs(U(2,1)), 0.0d0, tol, "observable: rectangle U column 1, y")
    call assert_equal_real(abs(U(3,1)), 0.0d0, tol, "observable: rectangle U column 1, z")
    call assert_equal_real(abs(U(1,2)), 0.0d0, tol, "observable: rectangle U column 2, x")
    call assert_equal_real(abs(U(2,2)), 1.0d0, tol, "observable: rectangle U column 2, y")
    call assert_equal_real(abs(U(3,2)), 0.0d0, tol, "observable: rectangle U column 2, z")
    call assert_equal_real(abs(U(1,3)), 0.0d0, tol, "observable: rectangle U column 3, x")
    call assert_equal_real(abs(U(2,3)), 0.0d0, tol, "observable: rectangle U column 3, y")
    call assert_equal_real(abs(U(3,3)), 1.0d0, tol, "observable: rectangle U column 3, z")
  end subroutine test_observable_full_rank_rectangle

  !> Three collinear points (intrinsic rank 1) embedded in a 5D ambient
  !| space: economy-mode rank = min(D,k) = 3 < D = 5, so U columns 4-5 and
  !| eigenvalues 4-5 must be observable's own zero-padding, not LAPACK
  !| output.
  subroutine test_observable_low_rank_padding()
    integer(int32), parameter :: d_dim = 5, n = 3
    real(real64)   :: vectors(d_dim, n)
    logical(c_bool)        :: member_selection_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr
    real(real64), parameter :: tol = 1.0d-8

    vectors      = 0.0d0
    vectors(1,1) = 0.0d0
    vectors(1,2) = 1.0d0
    vectors(1,3) = 2.0d0
    member_selection_mask = .true.

    call observable(vectors, d_dim, n, member_selection_mask, n, &
                    U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_int(d, 1, "observable: collinear-in-5D intrinsic dimension")
    call assert_equal_real(eigenvalues(1), 1.0d0, tol, "observable: collinear-in-5D lambda_1")
    call assert_equal_array_real(eigenvalues(4:5), [0.0d0, 0.0d0], 2, tol, &
                                 "observable: collinear-in-5D padded eigenvalues (4:5) must be exactly zero")
    call assert_equal_array_real(U(:,4), [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0], d_dim, tol, &
                                 "observable: collinear-in-5D padded U column 4 must be exactly zero")
    call assert_equal_array_real(U(:,5), [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0], d_dim, tol, &
                                 "observable: collinear-in-5D padded U column 5 must be exactly zero")
    call assert_equal_real(tangent_scales_value(1), 1.0d0, tol, "observable: collinear-in-5D tangent_scales(1)")
    call assert_equal_array_real(tangent_scales_value(2:5), [0.0d0,0.0d0,0.0d0,0.0d0], 4, tol, &
                                 "observable: collinear-in-5D padded tangent_scales(2:5) must be exactly zero")
    call assert_true(G > 1.0d10, "observable: collinear-in-5D spectral gap should be huge")
  end subroutine test_observable_low_rank_padding

  !> A single-member ensemble must be rejected: no meaningful covariance/SVD.
  !| Now caught automatically by the generated wrapper's n_selected_member
  !| DM_MIN(2) validation -- no bespoke kernel-level check needed.
  subroutine test_observable_too_few_members()
    integer(int32), parameter :: d_dim = 3, n = 4
    real(real64)   :: vectors(d_dim, n)
    logical(c_bool)        :: member_selection_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr

    vectors(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    vectors(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    vectors(:, 3) = [0.0d0, 1.0d0, 0.0d0]
    vectors(:, 4) = [2.0d0, 1.0d0, 0.0d0]
    member_selection_mask = .false.
    member_selection_mask(1) = .true.

    call observable(vectors, d_dim, n, member_selection_mask, 1_int32, &
                    U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "observable should reject an ensemble with fewer than 2 members")
  end subroutine test_observable_too_few_members

  !> n_dimensions=1 leaves no room for even a single spectral-gap
  !| comparison (r ranges 1..D-1) and must be rejected.
  subroutine test_observable_dimension_too_small()
    integer(int32), parameter :: d_dim = 1, n = 3
    real(real64)   :: vectors(d_dim, n)
    logical(c_bool)        :: member_selection_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr

    vectors(:, 1) = [0.0d0]
    vectors(:, 2) = [1.0d0]
    vectors(:, 3) = [2.0d0]
    member_selection_mask = .true.

    call observable(vectors, d_dim, n, member_selection_mask, n, &
                    U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "observable should reject n_dimensions=1")
  end subroutine test_observable_dimension_too_small

  ! --- ensemble_final_observable ------------------------------------------

  !> D=2, o=2, one ensemble: column 1 is the true accepted state (d=1, G=1.0,
  !| mu=[0.5,0.0]); column 2 is a rejected candidate with deliberately different
  !| values (d=0, G=99.0, mu=[9,9], accepted_history=.false.), mirroring
  !| `stc_push_ensemble_history`'s own trailing-rejected-push behavior (see
  !| `test_json_rejected_trailing_column_uses_last_accepted` for the same scenario at
  !| the JSON layer). The extraction must land on column 1, never column 2.
  subroutine test_ensemble_final_observable_trailing_rejected_column()
    integer(int32), parameter :: d_dim = 2, o = 2, n_e = 1
    real(real64)   :: U_hist(d_dim, d_dim, o, n_e), S_hist(d_dim, o, n_e), mu_hist(d_dim, o, n_e)
    real(real64)   :: G_hist(o, n_e)
    integer(int32) :: d_hist(o, n_e), k_hist(o, n_e)
    logical(c_bool)        :: accepted_hist(o, n_e)
    real(real64)   :: U_final(d_dim, d_dim, n_e), S_final(d_dim, n_e), mu_final(d_dim, n_e)
    real(real64)   :: G_final(n_e)
    integer(int32) :: d_final(n_e), k_final(n_e), final_index(n_e)
    logical(c_bool)        :: has_final(n_e)
    integer(int32) :: ierr

    k_hist(1, 1) = 2; k_hist(2, 1) = 3
    d_hist(1, 1) = 1; d_hist(2, 1) = 0
    G_hist(1, 1) = 1.0d0; G_hist(2, 1) = 99.0d0
    mu_hist(:, 1, 1) = [0.5d0, 0.0d0]; mu_hist(:, 2, 1) = [9.0d0, 9.0d0]
    S_hist(:, 1, 1) = [0.5d0, 2.0d0]; S_hist(:, 2, 1) = [7.0d0, 7.0d0]
    U_hist(:, 1, 1, 1) = [1.0d0, 0.0d0]; U_hist(:, 2, 1, 1) = [0.0d0, 1.0d0]
    U_hist(:, :, 2, 1) = 0.0d0
    accepted_hist(1, 1) = .true.
    accepted_hist(2, 1) = .false.

    call ensemble_final_observable(d_dim, o, n_e, U_hist, d_hist, S_hist, mu_hist, G_hist, k_hist, accepted_hist, &
                                   U_final, d_final, S_final, mu_final, G_final, k_final, has_final, final_index, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'ensemble_final_observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_true(has_final(1), "trailing rejected: has_final must be true")
    call assert_equal_int(final_index(1), 1, "trailing rejected: final_index must be column 1, not 2")
    call assert_equal_int(d_final(1), 1, "trailing rejected: d from column 1, not 0")
    call assert_equal_int(k_final(1), 2, "trailing rejected: k from column 1, not 3")
    call assert_equal_real(G_final(1), 1.0d0, 1.0d-12, "trailing rejected: G from column 1, not 99.0")
    call assert_equal_array_real(mu_final(:, 1), [0.5d0, 0.0d0], d_dim, 1.0d-12, &
                                 "trailing rejected: mu from column 1, not [9,9]")
    call assert_equal_array_real(S_final(:, 1), [0.5d0, 2.0d0], d_dim, 1.0d-12, &
                                 "trailing rejected: S from column 1, not [7,7]")
    call assert_equal_array_real(U_final(:, 1, 1), [1.0d0, 0.0d0], d_dim, 1.0d-12, &
                                 "trailing rejected: U column 1 from history column 1")
  end subroutine test_ensemble_final_observable_trailing_rejected_column

  !> Both history columns accepted (the common case): extraction must land on the
  !| last (most recent) populated column, column 2.
  subroutine test_ensemble_final_observable_no_rejection()
    integer(int32), parameter :: d_dim = 2, o = 2, n_e = 1
    real(real64)   :: U_hist(d_dim, d_dim, o, n_e), S_hist(d_dim, o, n_e), mu_hist(d_dim, o, n_e)
    real(real64)   :: G_hist(o, n_e)
    integer(int32) :: d_hist(o, n_e), k_hist(o, n_e)
    logical(c_bool)        :: accepted_hist(o, n_e)
    real(real64)   :: U_final(d_dim, d_dim, n_e), S_final(d_dim, n_e), mu_final(d_dim, n_e)
    real(real64)   :: G_final(n_e)
    integer(int32) :: d_final(n_e), k_final(n_e), final_index(n_e)
    logical(c_bool)        :: has_final(n_e)
    integer(int32) :: ierr

    k_hist(1, 1) = 2; k_hist(2, 1) = 3
    d_hist(1, 1) = 0; d_hist(2, 1) = 1
    G_hist(1, 1) = 2.0d0; G_hist(2, 1) = 1.5d0
    mu_hist(:, 1, 1) = [0.5d0, 0.0d0]; mu_hist(:, 2, 1) = [1.0d0, 0.0d0]
    S_hist(:, 1, 1) = [2.0d0, 0.0d0]; S_hist(:, 2, 1) = [0.5d0, 3.0d0]
    U_hist(:, 1, 1, 1) = [1.0d0, 0.0d0]; U_hist(:, 2, 1, 1) = [0.0d0, 1.0d0]
    U_hist(:, 1, 2, 1) = [1.0d0, 0.0d0]; U_hist(:, 2, 2, 1) = [0.0d0, 1.0d0]
    accepted_hist = .true.

    call ensemble_final_observable(d_dim, o, n_e, U_hist, d_hist, S_hist, mu_hist, G_hist, k_hist, accepted_hist, &
                                   U_final, d_final, S_final, mu_final, G_final, k_final, has_final, final_index, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'ensemble_final_observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_true(has_final(1), "no rejection: has_final must be true")
    call assert_equal_int(final_index(1), 2, "no rejection: final_index must be the last populated column")
    call assert_equal_int(d_final(1), 1, "no rejection: d from column 2")
    call assert_equal_int(k_final(1), 3, "no rejection: k from column 2")
    call assert_equal_real(G_final(1), 1.5d0, 1.0d-12, "no rejection: G from column 2")
    call assert_equal_array_real(mu_final(:, 1), [1.0d0, 0.0d0], d_dim, 1.0d-12, "no rejection: mu from column 2")
  end subroutine test_ensemble_final_observable_no_rejection

  !> Every history column has k=0 (unpopulated) -- no SVD ever ran for this ensemble
  !| (only possible for STOP_REASON_MAX_SIZE firing at the bootstrap step itself, see
  !| the kernel's own doc comment). has_final must be false and every _final output
  !| must come back exactly zero, not garbage.
  subroutine test_ensemble_final_observable_has_final_false()
    integer(int32), parameter :: d_dim = 2, o = 2, n_e = 1
    real(real64)   :: U_hist(d_dim, d_dim, o, n_e), S_hist(d_dim, o, n_e), mu_hist(d_dim, o, n_e)
    real(real64)   :: G_hist(o, n_e)
    integer(int32) :: d_hist(o, n_e), k_hist(o, n_e)
    logical(c_bool)        :: accepted_hist(o, n_e)
    real(real64)   :: U_final(d_dim, d_dim, n_e), S_final(d_dim, n_e), mu_final(d_dim, n_e)
    real(real64)   :: G_final(n_e)
    integer(int32) :: d_final(n_e), k_final(n_e), final_index(n_e)
    logical(c_bool)        :: has_final(n_e)
    integer(int32) :: ierr

    k_hist = 0
    d_hist = 0
    G_hist = 0.0d0
    mu_hist = 0.0d0
    S_hist = 0.0d0
    U_hist = 0.0d0
    accepted_hist = .false.

    call ensemble_final_observable(d_dim, o, n_e, U_hist, d_hist, S_hist, mu_hist, G_hist, k_hist, accepted_hist, &
                                   U_final, d_final, S_final, mu_final, G_final, k_final, has_final, final_index, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'ensemble_final_observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_true(.not. has_final(1), "no history: has_final must be false")
    call assert_equal_int(final_index(1), 0, "no history: final_index must be 0")
    call assert_equal_int(d_final(1), 0, "no history: d_final must be zero")
    call assert_equal_int(k_final(1), 0, "no history: k_final must be zero")
    call assert_equal_real(G_final(1), 0.0d0, 1.0d-12, "no history: G_final must be zero")
    call assert_equal_array_real(mu_final(:, 1), [0.0d0, 0.0d0], d_dim, 1.0d-12, "no history: mu_final must be zero")
    call assert_equal_array_real(S_final(:, 1), [0.0d0, 0.0d0], d_dim, 1.0d-12, "no history: S_final must be zero")
    call assert_equal_array_real(U_final(:, 1, 1), [0.0d0, 0.0d0], d_dim, 1.0d-12, "no history: U_final must be zero")
  end subroutine test_ensemble_final_observable_has_final_false

  !> Two ensembles with independent, deliberately-different histories in the same call --
  !| ensemble 1 has a trailing rejected column (must resolve to its column 1), ensemble 2
  !| has no rejection at all (must resolve to its last column, column 2). Proves the
  !| per-ensemble scan does not leak state across the ensemble dimension.
  subroutine test_ensemble_final_observable_multi_ensemble_independence()
    integer(int32), parameter :: d_dim = 2, o = 2, n_e = 2
    real(real64)   :: U_hist(d_dim, d_dim, o, n_e), S_hist(d_dim, o, n_e), mu_hist(d_dim, o, n_e)
    real(real64)   :: G_hist(o, n_e)
    integer(int32) :: d_hist(o, n_e), k_hist(o, n_e)
    logical(c_bool)        :: accepted_hist(o, n_e)
    real(real64)   :: U_final(d_dim, d_dim, n_e), S_final(d_dim, n_e), mu_final(d_dim, n_e)
    real(real64)   :: G_final(n_e)
    integer(int32) :: d_final(n_e), k_final(n_e), final_index(n_e)
    logical(c_bool)        :: has_final(n_e)
    integer(int32) :: ierr

    ! Ensemble 1: column 1 accepted, column 2 rejected.
    k_hist(1, 1) = 2; k_hist(2, 1) = 3
    d_hist(1, 1) = 1; d_hist(2, 1) = 0
    G_hist(1, 1) = 1.0d0; G_hist(2, 1) = 99.0d0
    mu_hist(:, 1, 1) = [0.5d0, 0.0d0]; mu_hist(:, 2, 1) = [9.0d0, 9.0d0]
    S_hist(:, 1, 1) = [0.5d0, 2.0d0]; S_hist(:, 2, 1) = [7.0d0, 7.0d0]
    U_hist(:, 1, 1, 1) = [1.0d0, 0.0d0]; U_hist(:, 2, 1, 1) = [0.0d0, 1.0d0]
    U_hist(:, :, 2, 1) = 0.0d0
    accepted_hist(1, 1) = .true.
    accepted_hist(2, 1) = .false.

    ! Ensemble 2: both columns accepted -- final must be column 2.
    k_hist(1, 2) = 4; k_hist(2, 2) = 5
    d_hist(1, 2) = 0; d_hist(2, 2) = 1
    G_hist(1, 2) = 3.0d0; G_hist(2, 2) = 2.5d0
    mu_hist(:, 1, 2) = [2.0d0, 2.0d0]; mu_hist(:, 2, 2) = [3.0d0, 3.0d0]
    S_hist(:, 1, 2) = [1.0d0, 0.0d0]; S_hist(:, 2, 2) = [1.0d0, 4.0d0]
    U_hist(:, 1, 1, 2) = [1.0d0, 0.0d0]; U_hist(:, 2, 1, 2) = [0.0d0, 1.0d0]
    U_hist(:, 1, 2, 2) = [0.0d0, 1.0d0]; U_hist(:, 2, 2, 2) = [1.0d0, 0.0d0]
    accepted_hist(:, 2) = .true.

    call ensemble_final_observable(d_dim, o, n_e, U_hist, d_hist, S_hist, mu_hist, G_hist, k_hist, accepted_hist, &
                                   U_final, d_final, S_final, mu_final, G_final, k_final, has_final, final_index, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'ensemble_final_observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_true(has_final(1), "multi-ensemble: ensemble 1 has_final must be true")
    call assert_equal_int(final_index(1), 1, "multi-ensemble: ensemble 1 resolves to its own column 1")
    call assert_equal_real(G_final(1), 1.0d0, 1.0d-12, "multi-ensemble: ensemble 1 G from its own column 1")

    call assert_true(has_final(2), "multi-ensemble: ensemble 2 has_final must be true")
    call assert_equal_int(final_index(2), 2, "multi-ensemble: ensemble 2 resolves to its own column 2")
    call assert_equal_real(G_final(2), 2.5d0, 1.0d-12, "multi-ensemble: ensemble 2 G from its own column 2")
    call assert_equal_int(k_final(2), 5, "multi-ensemble: ensemble 2 k from its own column 2")
    call assert_equal_array_real(mu_final(:, 2), [3.0d0, 3.0d0], d_dim, 1.0d-12, &
                                 "multi-ensemble: ensemble 2 mu from its own column 2, not ensemble 1's")
  end subroutine test_ensemble_final_observable_multi_ensemble_independence

  !> o=1: the single available column holds a rejected candidate. There is no earlier
  !| accepted column at all to fall back to (the window is too small to have retained
  !| one) -- has_final must be false, exactly as the kernel's own doc comment describes
  !| for "a small o lets a rejected push evict every accepted entry the window ever held".
  subroutine test_ensemble_final_observable_small_o_evicts_accepted()
    integer(int32), parameter :: d_dim = 2, o = 1, n_e = 1
    real(real64)   :: U_hist(d_dim, d_dim, o, n_e), S_hist(d_dim, o, n_e), mu_hist(d_dim, o, n_e)
    real(real64)   :: G_hist(o, n_e)
    integer(int32) :: d_hist(o, n_e), k_hist(o, n_e)
    logical(c_bool)        :: accepted_hist(o, n_e)
    real(real64)   :: U_final(d_dim, d_dim, n_e), S_final(d_dim, n_e), mu_final(d_dim, n_e)
    real(real64)   :: G_final(n_e)
    integer(int32) :: d_final(n_e), k_final(n_e), final_index(n_e)
    logical(c_bool)        :: has_final(n_e)
    integer(int32) :: ierr

    k_hist(1, 1) = 3
    d_hist(1, 1) = 0
    G_hist(1, 1) = 99.0d0
    mu_hist(:, 1, 1) = [9.0d0, 9.0d0]
    S_hist(:, 1, 1) = [7.0d0, 7.0d0]
    U_hist(:, :, 1, 1) = 0.0d0
    accepted_hist(1, 1) = .false.

    call ensemble_final_observable(d_dim, o, n_e, U_hist, d_hist, S_hist, mu_hist, G_hist, k_hist, accepted_hist, &
                                   U_final, d_final, S_final, mu_final, G_final, k_final, has_final, final_index, ierr)
    if (.not. is_ok(ierr)) then
      write (*, *) 'ensemble_final_observable failed unexpectedly: ', ierr
      error stop
    end if

    call assert_true(.not. has_final(1), "small o evicts accepted: has_final must be false")
    call assert_equal_int(final_index(1), 0, "small o evicts accepted: final_index must be 0")
    call assert_equal_real(G_final(1), 0.0d0, 1.0d-12, "small o evicts accepted: G_final must be zero, not 99.0")
  end subroutine test_ensemble_final_observable_small_o_evicts_accepted

end module mod_test_shape_truthful_clustering_observable
